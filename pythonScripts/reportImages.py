"""Insert/replace PHASE figures without duplicating images on report repair."""
from openpyxl.drawing.image import Image
from openpyxl.utils.cell import get_column_letter


def add_report_image(sheet, path, position, width, height):
    image = Image(path)  # Validate before replacing an existing figure.
    image.width, image.height = int(width), int(height)

    def cell(existing):
        anchor = existing.anchor
        if isinstance(anchor, str):
            return anchor.upper()
        marker = getattr(anchor, "_from", None)
        if marker is not None:
            return f"{get_column_letter(marker.col + 1)}{marker.row + 1}"
        return None

    sheet._images[:] = [existing for existing in sheet._images if cell(existing) != position.upper()]
    sheet.add_image(image, position)
