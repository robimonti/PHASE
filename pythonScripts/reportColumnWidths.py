"""Size report columns from displayed numbers, without changing their values."""
from datetime import date, datetime
from numbers import Real
import re
from math import ceil


def autofit_numeric_columns(sheet):
    # Report tables can be transposed: coordinates are not always columns D/E.
    # Ignore prose and figure captions, which would make the whole sheet wide.
    for cells in sheet.iter_cols():
        width = 0
        for cell in cells:
            value = cell.value
            if isinstance(value, bool) or value is None:
                continue
            if isinstance(value, Real):
                match = re.search(r"0\.([0#]+)", cell.number_format)
                if match:
                    text = f"{value:.{len(match[1])}f}"
                else:
                    text = str(value)
            elif isinstance(value, (datetime, date)):
                text = value.strftime("%d-%b-%Y")
            else:
                continue
            # Calibri 12 plus padding; include minus signs and decimal places.
            width = max(width, ceil(len(text) * 1.15 + 3))
        if width:
            column = sheet.column_dimensions[cells[0].column_letter]
            column.width = max(column.width or 0, width)
