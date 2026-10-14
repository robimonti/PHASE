import sys
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill

def format_excel_sheet5(report_path, sheet_name):
    # Load the workbook and select the sheet
    wb = openpyxl.load_workbook(report_path)
    sheet = wb[sheet_name]

    # General formatting for the table
    subtitle_row = 5  # Start of the table (below figure space)
    for row in sheet.iter_rows(min_row=subtitle_row+1, max_row=sheet.max_row, min_col=1, max_col=sheet.max_column):
        for cell in row:
            if cell.value is not None:
                cell.font = Font(size=12)
                cell.alignment = Alignment(horizontal='left', vertical='center')
    # Info
    for row in sheet.iter_rows(min_row=subtitle_row+2, max_row=sheet.max_row, min_col=2, max_col=5):
        for cell in row:
            if cell.value is not None:
                cell.alignment = Alignment(horizontal='center', vertical='center')
                cell.number_format = '0.000'  # Three decimals
    # Model-derived contrast, displacement, and persistence metrics.
    for row in sheet.iter_rows(min_row=subtitle_row+2, max_row=sheet.max_row, min_col=6, max_col=9):
        for cell in row:
            if cell.value is not None:
                cell.alignment = Alignment(horizontal='center', vertical='center')
                cell.number_format = '0.0'  # One decimals

    # The classification is descriptive, not a risk percentage.
    fills = {
        "typical": PatternFill(start_color="DCEFE5", end_color="DCEFE5", fill_type="solid"),
        "transient": PatternFill(start_color="FFF1CF", end_color="FFF1CF", fill_type="solid"),
        "persistent": PatternFill(start_color="FADBDD", end_color="FADBDD", fill_type="solid"),
    }
    for row in sheet.iter_rows(min_row=subtitle_row+2, max_row=sheet.max_row, min_col=10, max_col=10):
        for cell in row:
            cell.font = Font(size=12, bold=True)
            cell.alignment = Alignment(horizontal='center', vertical='center')
            if cell.value is not None:
                cell.fill = fills.get(str(cell.value).lower(), fills["typical"])

    # Format the title (A1)
    sheet['A1'].font = Font(size=20, bold=True)
    sheet['A1'].alignment = Alignment(horizontal='left', vertical='center')
    # Format the first subtitle (A2)
    sheet['A2'].font = Font(size=16, bold=False)
    sheet['A2'].alignment = Alignment(horizontal='left', vertical='center')
    # Format the second subtitle (A3)
    sheet['A3'].font = Font(size=16, bold=False)
    sheet['A3'].alignment = Alignment(horizontal='left', vertical='center')

    # Format the fourth subtitle (A5)
    sheet['A5'].font = Font(size=18, bold=True)
    sheet['A5'].alignment = Alignment(horizontal='left', vertical='center')

    # Save the workbook
    wb.save(report_path)
    print(f"Formatted sheet '{sheet_name}' and updated Excel file: {report_path}")

if __name__ == "__main__":
    # Parse arguments from MATLAB
    report_file = sys.argv[1]
    sheet_name = sys.argv[2]

    # Format the second sheet
    format_excel_sheet5(report_file, sheet_name)
