"""Large coordinates must remain readable in normal and transposed reports."""
import importlib.util
import pytest


def test_report_width_uses_formatted_values_without_changing_numbers(phase_root):
    openpyxl = pytest.importorskip("openpyxl")
    spec = importlib.util.spec_from_file_location(
        "report_widths", phase_root / "pythonScripts/reportColumnWidths.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    book = openpyxl.Workbook()
    sheet = book.active
    sheet["A1"] = "A very long report heading " * 20
    sheet["A33"] = "Y (UTM)"
    sheet["B33"] = 9531449.427630786
    sheet["B33"].number_format = "0.000"
    sheet["C33"] = -9531449.427630786
    sheet["C33"].number_format = "0.000000"
    sheet["D33"] = 0
    sheet["D33"].number_format = "0.000"
    original = tuple(sheet.values)
    module.autofit_numeric_columns(sheet)
    assert sheet.column_dimensions["B"].width >= 16
    assert sheet.column_dimensions["C"].width > sheet.column_dimensions["B"].width
    assert sheet.column_dimensions["D"].width >= 8
    assert sheet.column_dimensions["A"].width < 25
    assert tuple(sheet.values) == original


def test_every_report_sheet_applies_widths_before_save(phase_root):
    for number in range(1, 8):
        source = (phase_root / f"pythonScripts/reportSheet{number}.py").read_text()
        assert source.index("autofit_numeric_columns(sheet)") < source.index("wb.save(report_path)")
