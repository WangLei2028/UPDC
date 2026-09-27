from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_english_release_contains_no_cjk_documentation():
    files = [
        ROOT / "README.md",
        ROOT / "CHANGELOG.md",
        ROOT / "RELEASE_NOTES.md",
        *sorted((ROOT / "docs").glob("*.md")),
        *sorted((ROOT / ".github").rglob("*.md")),
        *sorted((ROOT / ".github").rglob("*.yml")),
    ]
    for path in files:
        text = path.read_text(encoding="utf-8")
        assert not any("\u4e00" <= ch <= "\u9fff" for ch in text), path


def test_complete_english_manuals_present():
    required = [
        "UPDC_COMPLETE_MANUAL_EN.md",
        "UPDC_ALGORITHM_MANUAL_EN.md",
        "UPDC_MANUAL_EN.md",
    ]
    for name in required:
        assert (ROOT / "docs" / name).is_file(), name
