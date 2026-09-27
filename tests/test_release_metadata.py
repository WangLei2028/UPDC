from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_standard_mit_license_present():
    text = (ROOT / "LICENSE").read_text(encoding="utf-8")
    assert text.startswith("MIT License\n")
    assert "Permission is hereby granted, free of charge" in text
    assert "THE SOFTWARE IS PROVIDED \"AS IS\"" in text


def test_citation_metadata_contains_both_required_dois():
    cff_text = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert "doi: 10.1093/mnras/stac2664" in cff_text
    assert "doi: 10.3847/1538-4365/ad9566" in cff_text

    text = (ROOT / "CITATION.md").read_text(encoding="utf-8")
    assert "10.1093/mnras/stac2664" in text
    assert "10.3847/1538-4365/ad9566" in text


def test_no_prepublication_license_placeholder_remains():
    assert not (ROOT / "docs" / "LICENSE_BEFORE_PUBLICATION.md").exists()
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "grants no general open-source license" not in readme
