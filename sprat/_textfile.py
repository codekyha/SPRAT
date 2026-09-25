"""Parser and schema for the two plain-text input files (``.phc`` and ``.par``).

Grammar
-------
* ``[section]`` starts a section; section names are case-insensitive.
* ``key = value`` assigns a value; keys are case-insensitive, values keep their case.
* ``#`` starts a comment (to the end of the line); blank lines are ignored.
* A line without ``=`` inside a section that allows bare lines (``[rods]``) is kept verbatim.
* Values are typed by a schema.  Unknown sections or keys are errors, so a misspelled key can
  never be silently ignored.

Only the standard library is used.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Any


class TextFileError(ValueError):
    """Raised for any syntax or schema error in a ``.phc`` or ``.par`` file."""


@dataclass
class Field:
    """One key of a section: its type, default, documentation and constraints."""

    type: str                       # float | int | bool | str | enum | float_or_auto | float_or_calibrate | int_or_auto
    default: Any
    help: str
    choices: tuple = ()
    unit: str = ""
    minimum: float | None = None
    maximum: float | None = None
    required: bool = False


@dataclass
class Section:
    """A section of the schema.

    ``fields`` maps key -> Field.  ``free_keys`` allows arbitrary keys (typed by ``free_type``);
    ``bare_lines`` allows lines without ``=``.
    """

    help: str
    fields: dict = field(default_factory=dict)
    free_keys: bool = False
    free_type: str = "float"
    bare_lines: bool = False


# --------------------------------------------------------------------------- lexing
_SECTION_RE = re.compile(r"^\[\s*([A-Za-z0-9_.\-]+)\s*\]\s*$")


def strip_comment(line: str) -> str:
    i = line.find("#")
    return line if i < 0 else line[:i]


def parse_text(text: str, path: str = "<text>") -> dict[str, list[tuple[int, str | None, str]]]:
    """Return {section: [(line_number, key or None, value), ...]} preserving order.

    ``key`` is None for a bare line.  Keys are lower-cased.
    """
    sections: dict[str, list[tuple[int, str | None, str]]] = {}
    current: str | None = None
    for number, raw in enumerate(text.splitlines(), start=1):
        line = strip_comment(raw).strip()
        if not line:
            continue
        m = _SECTION_RE.match(line)
        if m:
            current = m.group(1).lower()
            sections.setdefault(current, [])
            continue
        if current is None:
            raise TextFileError(f"{path}:{number}: a value appears before any [section] header: {raw.strip()!r}")
        if "=" in line:
            key, _, value = line.partition("=")
            key = key.strip().lower()
            if not key:
                raise TextFileError(f"{path}:{number}: empty key")
            sections[current].append((number, key, value.strip()))
        else:
            sections[current].append((number, None, line))
    return sections


# --------------------------------------------------------------------------- typing
_TRUE = {"true", "yes", "on", "1"}
_FALSE = {"false", "no", "off", "0"}


def parse_scalar(text: str, ftype: str, where: str, fld: Field | None = None) -> Any:
    """Convert one textual value to the type of the field."""
    t = text.strip()
    if ftype == "float":
        try:
            v = float(t)
        except ValueError:
            raise TextFileError(f"{where}: expected a number, got {t!r}") from None
        if not math.isfinite(v):
            raise TextFileError(f"{where}: number must be finite, got {t!r}")
    elif ftype == "int":
        try:
            v = int(t)
        except ValueError:
            try:
                fv = float(t)
            except ValueError:
                raise TextFileError(f"{where}: expected an integer, got {t!r}") from None
            if fv != int(fv):
                raise TextFileError(f"{where}: expected an integer, got {t!r}")
            v = int(fv)
    elif ftype == "bool":
        low = t.lower()
        if low in _TRUE:
            v = True
        elif low in _FALSE:
            v = False
        else:
            raise TextFileError(f"{where}: expected true/false, got {t!r}")
    elif ftype == "str":
        v = t
    elif ftype == "enum":
        low = t.lower()
        choices = tuple(fld.choices) if fld else ()
        if low not in choices:
            raise TextFileError(f"{where}: expected one of {', '.join(choices)}, got {t!r}")
        v = low
    elif ftype == "float_or_auto":
        v = "auto" if t.lower() == "auto" else parse_scalar(t, "float", where)
    elif ftype == "int_or_auto":
        v = "auto" if t.lower() == "auto" else parse_scalar(t, "int", where)
    elif ftype == "float_or_calibrate":
        v = "calibrate" if t.lower() == "calibrate" else parse_scalar(t, "float", where)
    else:
        raise TextFileError(f"{where}: unknown field type {ftype!r} (schema bug)")
    if fld is not None and isinstance(v, (int, float)) and not isinstance(v, bool):
        if fld.minimum is not None and v < fld.minimum:
            raise TextFileError(f"{where}: {v} is below the minimum {fld.minimum}")
        if fld.maximum is not None and v > fld.maximum:
            raise TextFileError(f"{where}: {v} is above the maximum {fld.maximum}")
    return v


def parse_list(text: str, ftype: str, where: str, fld: Field | None = None) -> list:
    """``a:b:step`` (inclusive range) or ``v1, v2, v3`` -> list of typed values."""
    t = text.strip()
    if ":" in t and "," not in t:
        parts = [p.strip() for p in t.split(":")]
        if len(parts) != 3:
            raise TextFileError(f"{where}: a range needs start:stop:step, got {t!r}")
        a, b, h = (float(p) for p in parts)
        if h == 0 or (b - a) * h < 0:
            raise TextFileError(f"{where}: range {t!r} has an inconsistent step")
        n = int(round((b - a) / h)) + 1
        if abs(a + (n - 1) * h - b) > 1e-9 * max(1.0, abs(b)):
            raise TextFileError(f"{where}: range {t!r}: stop is not start + k*step")
        values = [a + i * h for i in range(n)]
        out = []
        for v in values:
            # round to the precision of the step so that 0.05 + 3*0.005 prints as 0.065, not 0.06500000001
            decimals = max(0, -int(math.floor(math.log10(abs(h)))) + 3) if h != 0 else 6
            vr = round(v, decimals)
            out.append(parse_scalar(repr(vr), ftype, where, fld))
        return out
    return [parse_scalar(p, ftype, where, fld) for p in t.split(",") if p.strip()]


# --------------------------------------------------------------------------- schema application
def apply_schema(sections: dict, schema: dict[str, Section], path: str) -> dict[str, Any]:
    """Validate the parsed sections against the schema and return a typed nested dict.

    Missing keys take their defaults; unknown sections or keys raise TextFileError with the
    list of valid names; required keys that are missing raise as well.
    """
    out: dict[str, Any] = {}
    for name in sections:
        if name not in schema:
            raise TextFileError(f"{path}: unknown section [{name}]; valid sections: "
                                + ", ".join(f"[{s}]" for s in schema))
    for name, sec in schema.items():
        entries = sections.get(name, [])
        values: dict[str, Any] = {k: f.default for k, f in sec.fields.items()}
        bare: list[str] = []
        seen: set[str] = set()
        for number, key, value in entries:
            where = f"{path}:{number} [{name}]"
            if key is None:
                if not sec.bare_lines:
                    raise TextFileError(f"{where}: expected key = value, got {value!r}")
                bare.append(value)
                continue
            if key in seen:
                raise TextFileError(f"{where}: key {key!r} given twice")
            seen.add(key)
            if key in sec.fields:
                fld = sec.fields[key]
                values[key] = parse_scalar(value, fld.type, f"{where} {key}", fld)
            elif sec.free_keys:
                values[key] = parse_scalar(value, sec.free_type, f"{where} {key}")
            else:
                raise TextFileError(f"{where}: unknown key {key!r}; valid keys: "
                                    + ", ".join(sec.fields) if sec.fields else f"{where}: section takes no keys")
        for key, fld in sec.fields.items():
            if fld.required and key not in seen:
                raise TextFileError(f"{path}: [{name}] {key} is required")
        if sec.bare_lines:
            values["_lines"] = bare
        out[name] = values
    return out


def schema_markdown(schema: dict[str, Section], title: str) -> str:
    """Render the schema as a Markdown reference (used to build the docs pages)."""
    lines = [f"# {title}", ""]
    for name, sec in schema.items():
        lines += [f"## [{name}]", "", sec.help, ""]
        if sec.fields:
            lines += ["| key | type | default | unit | meaning |", "|---|---|---|---|---|"]
            for key, f in sec.fields.items():
                typ = f.type if f.type != "enum" else " / ".join(f.choices)
                dflt = "required" if f.required else repr(f.default)
                lines.append(f"| `{key}` | {typ} | {dflt} | {f.unit} | {f.help} |")
            lines.append("")
        if sec.free_keys:
            lines += [f"Free keys allowed (each value is a {sec.free_type}).", ""]
        if sec.bare_lines:
            lines += ["Bare lines (without `=`) allowed.", ""]
    return "\n".join(lines)


def dotted_get(d: dict, dotted: str) -> Any:
    cur: Any = d
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            raise KeyError(dotted)
        cur = cur[part]
    return cur


def dotted_set(d: dict, dotted: str, value: Any) -> None:
    parts = dotted.split(".")
    cur = d
    for part in parts[:-1]:
        cur = cur.setdefault(part, {})
    cur[parts[-1]] = value
