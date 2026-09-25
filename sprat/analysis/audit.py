"""``sprat audit``: recompute the registry from the records and compare it, row by row, with an expected registry
(the one deposited with the paper, or the manuscript registry with its later sections).

Rows.  Every section of the expected registry is read:

* value rows (``quantity`` and ``value``) are record-derived; each is compared with the recomputed row of the same
  quantity name, wherever SPRAT's registry holds it.  Measured rows are also compared on the record they name and on
  the record fields they carry (f_r, margin, raw wavelength);
* superseded rows (``v3_value``, ``v4_value``, ``reason``) state which earlier value a current one replaced; they are
  matched to the recomputed statement by the text of ``v3_value`` and every number in their text is compared.  They
  are tallied apart from the value rows (verdict CONSISTENT);
* literature values, references and notes (rows with a ``key``, citation, DOI or URL, plain strings, and top-level
  text such as ``dataset_note``) are not record-derived: they are listed and counted separately, never as OK.

Tolerance: the printed precision of the registry value.  A number agrees when the recomputed value lies within one
unit of the last digit the registry prints: 1.9595 within 0.0001, 5663.6 within 0.1, "9.1e3" within 100, "34 %"
within 1.  Trailing zeros of a number printed without a decimal fraction are placeholders (97 500 was rounded to three
significant figures, and JSON writes it 97500.0), so 97500 is compared within 100 and 9000 within 1000.  Counts (the
dataset section) must agree exactly and text must agree exactly.  A recomputed value that agrees but would round to a
different last digit is marked "last digit".  Every deviation is listed; nothing is tuned.
"""

from __future__ import annotations

import json
import os
import re
import time
from typing import Any

from . import pipeline

METADATA = {"spec", "created", "rule", "software", "provenance"}
NUM = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
_FULL = re.compile(r"\s*([+-]?)(\d+)(?:\.(\d+))?(?:[eE]([+-]?\d+))?\s*(%|per cent)?\s*")
SLACK = 1e-9                      # relative slack of the comparison (floating-point representation of one unit)
RECORD_FIELDS = (("f_r", "f_r"), ("margin", "margin"), ("lambda_r_nm_raw", "lambda_r_nm_raw"))
IDENTITY_FIELDS = (("etiket", "label"), ("label", "label"), ("is_no", "job_id"), ("kod_sha256", "code_sha256"))


# --------------------------------------------------------------------------- numbers and their printed precision
def _text_of(v) -> str | None:
    if isinstance(v, bool) or v is None:
        return None
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        s = repr(v)
        return s[:-2] if s.endswith(".0") else s
    if isinstance(v, str):
        return v
    return None


def number(v) -> tuple[float, float] | None:
    """(value, printed unit) of a registry number or of a numeric string such as '9.1e3', '34 %' or '0.003 %'."""
    t = _text_of(v)
    if t is None:
        return None
    m = _FULL.fullmatch(t)
    if not m:
        return None
    sign, ip, frac, exp = m.group(1), m.group(2), m.group(3), int(m.group(4) or 0)
    val = float(f"{sign}{ip}{'.' + frac if frac else ''}e{exp}")
    if frac:
        unit = 10.0 ** (exp - len(frac))
    else:
        z = len(ip) - len(ip.rstrip("0")) if int(ip) != 0 else 0
        unit = 10.0 ** (exp + z)
    return val, unit


def _norm(s: str) -> str:
    return " ".join(str(s).split())


def compare(a, b, exact: bool = False) -> tuple[bool, str, bool]:
    """Compare a recomputed value ``a`` with a registry value ``b``: (agrees, detail, last digit differs)."""
    pb = number(b)
    if pb is not None and not isinstance(b, bool):
        pa = number(a) if not isinstance(a, (int, float)) or isinstance(a, bool) else (float(a), 0.0)
        if pa is None:
            return False, f"{a!r} vs {b!r}", False
        x, (y, u) = pa[0], pb
        if exact:
            return x == y, f"{a} vs {b} (exact)", False
        d = abs(x - y)
        ok = d <= u * (1 + SLACK) + 1e-300
        last = ok and d > 0.5 * u * (1 + SLACK)
        return ok, f"{a} vs {b} (tol {u:g})" + (" [last digit]" if last else ""), last
    if isinstance(b, list):
        if not isinstance(a, list) or len(a) != len(b):
            return False, f"{a!r} vs {b!r}", False
        r = [compare(x, y, exact) for x, y in zip(a, b)]
        return all(o for o, _, _ in r), "; ".join(m for _, m, _ in r), any(t for _, _, t in r)
    if isinstance(b, dict):
        if not isinstance(a, dict) or set(a) != set(b):
            return False, f"{a!r} vs {b!r}", False
        r = [compare(a[k], b[k], exact) for k in b]
        return all(o for o, _, _ in r), "; ".join(m for _, m, _ in r), any(t for _, _, t in r)
    return _norm(a) == _norm(b), f"{a!r} vs {b!r}", False


def _skeleton(s: str) -> str:
    return _norm(NUM.sub("#", s))


def compare_text(a: str, b: str) -> tuple[bool, str, bool]:
    """Two statements agree when their text around the numbers is identical and every number agrees."""
    if _skeleton(a) != _skeleton(b):
        return False, f"text differs: {a!r} vs {b!r}", False
    r = [compare(x, y) for x, y in zip(NUM.findall(a), NUM.findall(b))]
    return all(o for o, _, _ in r), "; ".join(m for _, m, _ in r), any(t for _, _, t in r)


# --------------------------------------------------------------------------- classification
def _kind_of(sec: str, row) -> str:
    if isinstance(row, str):
        return "reference" if sec.startswith("references") else "note"
    if isinstance(row, dict):
        if "quantity" in row and "value" in row and row.get("kind") not in ("literature", "reference"):
            return "value"
        if "v3_value" in row or "v4_value" in row:
            return "superseded"
        if "values" in row or sec.startswith("literature") or row.get("kind") == "literature":
            return "literature"
        if any(k in row for k in ("doi", "citation", "url", "key")) or sec.startswith("references"):
            return "reference"
    return "note"


def _index(new: dict) -> dict[str, tuple[str, dict]]:
    idx: dict[str, tuple[str, dict]] = {}
    for sec, rows in new.items():
        if isinstance(rows, list):
            for r in rows:
                if isinstance(r, dict) and "quantity" in r:
                    idx.setdefault(r["quantity"], (sec, r))
    return idx


def compare_registries(new: dict, expected: dict, rel: float | None = None) -> dict[str, Any]:
    """Audit every row of ``expected`` against the recomputed registry ``new``.  ``rel`` is accepted for compatibility
    with earlier versions and ignored: the tolerance is the printed precision of each registry value."""
    idx = _index(new)
    same_sec = {sec: {r["quantity"]: r for r in rows if isinstance(r, dict) and "quantity" in r}
                for sec, rows in new.items() if isinstance(rows, list)}
    sup_new = [r for r in (new.get("superseded") or []) if isinstance(r, dict) and "v3_value" in r]
    rows, sup, notrec = [], [], []
    for sec, items in expected.items():
        if sec in METADATA:
            continue
        if isinstance(items, str):
            notrec.append(dict(section=sec, item=sec, kind="note", text=items))
            continue
        if not isinstance(items, list):
            continue
        for i, r in enumerate(items):
            kind = _kind_of(sec, r)
            if kind == "value":
                q = r["quantity"]
                n = same_sec.get(sec, {}).get(q) or (idx.get(q) or (None, None))[1]
                if n is None:
                    rows.append(dict(section=sec, quantity=q, expected=r.get("value"), recomputed=None, verdict="MISSING", detail="", last_digit=False))
                    continue
                ok, detail, last = compare(n.get("value"), r.get("value"), exact=(sec == "dataset"))
                if r.get("kind") == "measured":
                    for ek, nk in IDENTITY_FIELDS:
                        if r.get(ek) not in (None, "") and n.get(nk) not in (None, "") and str(r[ek]) != str(n[nk]):
                            ok, detail = False, detail + f"; {ek} {r[ek]!r} vs {n[nk]!r}"
                    for ek, nk in RECORD_FIELDS:
                        if r.get(ek) is not None and n.get(nk) is not None:
                            o2, d2, l2 = compare(n[nk], r[ek])
                            ok, detail, last = ok and o2, detail + f"; {ek} {d2}", last or l2
                rows.append(dict(section=sec, quantity=q, expected=r.get("value"), recomputed=n.get("value"), verdict="OK" if ok else "DEVIATION",
                                 detail=detail, last_digit=last))
            elif kind == "superseded":
                hit = [s for s in sup_new if _skeleton(s.get("v3_value", "")) == _skeleton(r.get("v3_value", ""))]
                if not hit:
                    sup.append(dict(section=sec, index=i, v3_value=r.get("v3_value"), v4_value=r.get("v4_value"), recomputed=None,
                                    verdict="MISSING", detail="", last_digit=False))
                    continue
                s = hit[0]
                res = [compare_text(s.get(k, ""), r.get(k, "")) for k in ("v3_value", "v4_value", "reason")]
                ok = all(o for o, _, _ in res)
                sup.append(dict(section=sec, index=i, v3_value=r.get("v3_value"), v4_value=r.get("v4_value"),
                                recomputed=dict(v3_value=s.get("v3_value"), v4_value=s.get("v4_value"), reason=s.get("reason")),
                                verdict="CONSISTENT" if ok else "DEVIATION", detail=" | ".join(m for _, m, _ in res if m),
                                last_digit=any(t for _, _, t in res)))
            else:
                label = r if isinstance(r, str) else (r.get("key") or r.get("quantity") or json.dumps(r)[:60])
                notrec.append(dict(section=sec, item=label, kind=kind))
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    sc: dict[str, int] = {}
    for r in sup:
        sc[r["verdict"]] = sc.get(r["verdict"], 0) + 1
    nk: dict[str, int] = {}
    for r in notrec:
        nk[r["kind"]] = nk.get(r["kind"], 0) + 1
    return dict(generated=time.strftime("%Y-%m-%d %H:%M:%S"),
                tolerance="printed precision of each registry value (one unit of its last printed digit); counts and text exactly",
                rows=rows, counts=counts, superseded=sup, superseded_counts=sc, not_record_derived=notrec,
                not_record_derived_count=len(notrec), not_record_derived_kinds=nk,
                ok=not any(r["verdict"] != "OK" for r in rows) and not any(r["verdict"] != "CONSISTENT" for r in sup))


def _verdict_line(counts: dict[str, int], first: str) -> str:
    keys = sorted(counts, key=lambda k: (k != first, k))
    return ", ".join(f"{k} {counts[k]}" for k in keys) if counts else "none"


def markdown(cmp: dict) -> str:
    L = ["# Audit: recomputed registry against the expected one", "",
         f"generated {cmp['generated']}; tolerance: {cmp.get('tolerance', '')}", "",
         "| section | quantity | expected | recomputed | verdict |", "|---|---|---|---|---|"]
    for r in cmp["rows"]:
        v = r["verdict"] + (" (last digit)" if r.get("last_digit") else "")
        L.append(f"| {r['section']} | {r['quantity']} | {r['expected']} | {r['recomputed']} | **{v}** |")
    if cmp.get("superseded"):
        L += ["", "| section | earlier value (v3) | current value (v4) | verdict |", "|---|---|---|---|"]
        for r in cmp["superseded"]:
            L.append(f"| {r['section']} | {r['v3_value']} | {r['v4_value']} | **{r['verdict']}** |")
    if cmp.get("not_record_derived"):
        L += ["", "Not record-derived (listed, not audited):", ""]
        L += [f"- {r['section']}: {r['item'] if r['kind'] != 'note' else 'note'} ({r['kind']})" for r in cmp["not_record_derived"]]
    L += ["", "Verdicts: " + _verdict_line(cmp["counts"], "OK")]
    if cmp.get("superseded_counts"):
        L.append("Superseded statements: " + _verdict_line(cmp["superseded_counts"], "CONSISTENT"))
    if cmp.get("not_record_derived_count"):
        L.append("Not record-derived: %d (%s)" % (cmp["not_record_derived_count"],
                                                 ", ".join(f"{k} {v}" for k, v in sorted(cmp["not_record_derived_kinds"].items()))))
    L.append("")
    last = [r for r in cmp["rows"] + cmp.get("superseded", []) if r.get("last_digit")]
    if last:
        L += ["## Agreement in the last printed digit", "",
              "The recomputed value lies within one unit of the last digit but would round to a neighbouring value:", ""]
        L += [f"- {r.get('quantity') or r.get('v3_value')}: {r['detail']}" for r in last] + [""]
    dev = [r for r in cmp["rows"] + cmp.get("superseded", []) if r["verdict"] in ("DEVIATION", "MISSING")]
    if dev:
        L += ["## Deviations", ""] + [f"- {r.get('quantity') or r.get('v3_value')}: {r['verdict']} {r['detail']}" for r in dev] + [""]
    return "\n".join(L)


def audit_cli(records_dir: str, out: str, expected: str | None = None) -> None:
    A = pipeline.analyze(records_dir, out)
    if expected:
        with open(expected, encoding="utf-8") as fh:
            exp = json.load(fh)
        cmp = compare_registries(A["registry"], exp)
        cmp["expected_file"] = pipeline.portable_path(expected)
        with open(os.path.join(out, "audit.json"), "w", encoding="utf-8") as fh:
            json.dump(cmp, fh, indent=1)
        with open(os.path.join(out, "AUDIT_REPORT.md"), "w", encoding="utf-8") as fh:
            fh.write(markdown(cmp))
        print(markdown(cmp))
    else:
        print("no expected registry given; the recomputed registry is in", os.path.join(out, "numbers_registry.md"))


__all__ = ["number", "compare", "compare_text", "compare_registries", "markdown", "audit_cli"]
