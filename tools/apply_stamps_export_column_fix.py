"""Remove the artificial FZ column limit from the StaMPS MLAPP export."""

from pathlib import Path
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]
MLAPP = ROOT / "PHASE_Preprocessing" / "PHASE_StaMPS.mlapp"
OLD = "'Range', 'E1:FZ500000'"
OLD_WRAPPED = "'Range', 'D1:FZ500000'"


def main() -> None:
    with zipfile.ZipFile(MLAPP, "r") as source:
        entries = {name: source.read(name) for name in source.namelist()}
    document = entries["matlab/document.xml"].decode("utf-8")
    count = document.count(OLD) + document.count(OLD_WRAPPED)
    document = document.replace(OLD, "'Range', 'E1'")
    document = document.replace(OLD_WRAPPED, "'Range', 'D1'")
    entries["matlab/document.xml"] = document.encode("utf-8")

    with tempfile.NamedTemporaryFile(dir=MLAPP.parent, suffix=".mlapp", delete=False) as tmp:
        temporary = Path(tmp.name)
    try:
        with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as target:
            for name, data in entries.items():
                target.writestr(name, data)
        temporary.replace(MLAPP)
    finally:
        if temporary.exists():
            temporary.unlink()
    print(f"Removed {count} hard-coded StaMPS export column limits from {MLAPP}.")


if __name__ == "__main__":
    main()
