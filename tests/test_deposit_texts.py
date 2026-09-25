"""The texts that `sprat deposit` writes: the data availability statement, the metadata with the absorber spectra and
the preprint, the licence, and the version of the plane-wave record taken from the record itself. No word
"campaign" in any of them."""
from sprat import deposit

DOI = "10.5281/zenodo.22912910"
COUNTS = dict(harminv=858, valid=741, excluded_margin=17, excluded_gap=100, waveguide_mode_admitted=537, dead_zone_admitted=436,
              resolved=841, spectra=14, references=14, spectra_absorber=7, references_absorber=7, fields=2, bands=4, pwe=2)


def test_statement_and_counts():
    text = deposit.data_availability(COUNTS, data_doi=DOI, software_doi=DOI, raw_included=True, supplementary=True)
    assert "the scripts that produced the records" in text and "SPRAT version" in text and "campaign" not in text.lower()
    assert "| of these, with the guide continued into the absorber | 7 / 7 |" in text
    no_absorber = deposit.data_availability(dict(COUNTS, spectra=7, references=7, spectra_absorber=0, references_absorber=0))
    assert "absorber" not in no_absorber and "| spectra / references | 7 / 7 |" in no_absorber


def test_metadata_with_the_absorber_spectra_and_the_preprint():
    raw = dict(included=True, version="2.2.0", readme=True, changes="CHANGES_FROM_2.2.0.md")
    m = deposit.zenodo_metadata(COUNTS, "3.1.0", raw_included=True, included=["supplementary.pdf", "sprat-1.2.0.zip"], data_doi=DOI,
                                software_doi=DOI, org=dict(imported=893, sprat=0, sprat_pwe=1, sprat_pwe_version="1.1.0"), raw=raw,
                                manuscript_registry=True, h14b=True, arxiv_id="2609.12345")
    assert m["related_identifiers"] == [dict(relation="isSupplementTo", identifier="arXiv:2609.12345", scheme="arxiv",
                                             resource_type="publication-preprint")]
    d = m["description"]
    assert "7 with the guide ending in the perfectly matched layer, 7 with the guide continued into an absorber" in d
    assert "raw_h14b.tar.gz holds" in d and "arXiv:2609.12345" in d and "plane-wave layer of SPRAT 1.1.0" in d
    assert "apart from the changes listed in its CHANGES_FROM_2.2.0.md" in d
    assert "campaign" not in (d + m["notes"]).lower()
    mit = [x for x in m["licenses"] if x["id"] == "MIT"][0]
    assert mit["applies_to"] == ("the software SPRAT (sprat-1.2.0.zip), the scripts inside raw_legacy.tar.gz and the scripts inside "
                                 "raw_h14b.tar.gz")
    assert mit["applies_to"] + " under the MIT licence" in m["notes"]
    plain = deposit.zenodo_metadata(COUNTS, "3.1.0", data_doi=DOI, software_doi=DOI)
    assert "related_identifiers" not in plain


def test_licence_names_both_script_archives():
    text = deposit.licence_text(dict(included=True), dict(imported=1, sprat=0), sprat_archive="sprat-1.2.0.zip", h14b=True)
    assert "raw_legacy.tar.gz" in text and "raw_h14b.tar.gz" in text and "campaign" not in text.lower()


def test_plane_wave_record_keeps_its_version():
    recs = [dict(task=dict(source="legacy-campaign"), params=dict(run=dict(mode="harminv")), software=dict(version="1.2.0")),
            dict(task=dict(source="sprat"), params=dict(run=dict(mode="pwe")), software=dict(version="1.1.0"))]
    org = deposit.origin(recs)
    assert org == dict(imported=1, sprat=0, sprat_pwe=1, sprat_pwe_version="1.1.0")
    assert "plane-wave layer of SPRAT 1.1.0" in deposit.origin_text(org, raw_included=True, h14b=True)
    assert "`raw_legacy.tar.gz` and `raw_h14b.tar.gz`" in deposit.origin_text(org, raw_included=True, h14b=True)
