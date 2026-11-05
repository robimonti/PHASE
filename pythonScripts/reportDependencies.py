"""Check the actual report interpreter, including image encoding and XLSX saving."""
import argparse
import importlib
import io
from pathlib import Path
import subprocess
import sys


def missing_packages():
    missing = []
    for package, module in (("Pillow", "PIL.Image"), ("openpyxl", "openpyxl")):
        try:
            importlib.import_module(module)
        except ImportError:
            missing.append(package)
    return missing


def check_images():
    from PIL import Image
    from openpyxl import Workbook
    from openpyxl.drawing.image import Image as ExcelImage
    png = io.BytesIO()
    Image.new("RGB", (2, 2), "white").save(png, format="PNG")
    png.seek(0)
    book = Workbook()
    book.active.add_image(ExcelImage(png), "A1")
    book.save(io.BytesIO())
    book.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ensure", action="store_true", help="Install missing report packages in this interpreter")
    args = parser.parse_args()
    missing = missing_packages()
    if missing:
        if not args.ensure:
            raise RuntimeError(f"Missing report packages in {sys.executable}: {', '.join(missing)}")
        print(f"Preparing Excel report dependencies in {sys.executable}: {', '.join(missing)}", flush=True)
        subprocess.run([sys.executable, "-m", "pip", "install", "--disable-pip-version-check", "--no-input", *missing], check=True)
        # openpyxl may have cached an unavailable Pillow import. Validate in a
        # fresh process, just like the seven external report formatters.
        subprocess.run([sys.executable, str(Path(__file__).resolve())], check=True)
        return
    check_images()
    print(f"Excel report ready: openpyxl and Pillow image embedding verified ({sys.executable}).")


if __name__ == "__main__":
    main()
