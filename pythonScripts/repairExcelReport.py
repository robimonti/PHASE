"""Restore PHASE report figures from saved PNGs, without rerunning modelling."""
import argparse
from datetime import datetime
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def repair_report(report_path):
    import openpyxl
    from reportSheet1 import format_excel
    from reportSheet2 import format_excel_sheet2
    from reportSheet3 import format_excel_sheet3
    from reportSheet4 import format_excel_sheet4
    from reportSheet5 import format_excel_sheet5
    from reportSheet6 import format_excel_sheet6
    from reportSheet7 import format_excel_sheet7

    report = Path(report_path).resolve()
    figures = report.parent / 'figures'
    logo = Path(__file__).resolve().parents[1] / 'PHASE_logo.png'
    book = openpyxl.load_workbook(report)
    jobs = []
    expected = []
    if 'General' in book.sheetnames:
        images = [(str(logo), 'F1', 366, 180),
                  (str(figures / 'AOI_context.png'), 'E8', 900, 450),
                  (str(figures / 'AOI_PS.png'), 'E30', 900, 450)]
        expected.extend(Path(item[0]) for item in images)
        jobs.append((format_excel, [images, 'General']))
    for sheet, name, formatter in [
        ('Raw displacement', 'AOI_meanDispl_RAW.png', format_excel_sheet2),
        ('Modelled displacement', 'AOI_AvgVelocity.png', format_excel_sheet3),
        ('Uncertainty', 'AOI_uncertainty.png', format_excel_sheet4),
    ]:
        if sheet in book.sheetnames:
            path = figures / name
            expected.append(path)
            jobs.append((formatter, [str(path), sheet]))
    if 'AOI anomalies' in book.sheetnames:
        jobs.append((format_excel_sheet5, ['AOI anomalies']))
    if 'Interpolation' in book.sheetnames:
        sheet = book['Interpolation']
        paths = []
        for header in range(6, sheet.max_row + 1, 27):
            if sheet.cell(header, 1).value != 'ID':
                raise ValueError(f'Unrecognized query-point block at row {header}; workbook left unchanged.')
            identifier = sheet.cell(header + 1, 1).value
            path = figures / f'ExtrDispl_{identifier}.png'
            paths.append(str(path))
            expected.append(path)
        jobs.append((format_excel_sheet6, ['Interpolation', paths]))
    if 'Modelled PS' in book.sheetnames:
        jobs.append((format_excel_sheet7, ['Modelled PS']))
    book.close()
    if not jobs:
        raise ValueError('No recognized PHASE report sheets; workbook left unchanged.')
    missing = [str(path) for path in expected if not path.is_file()]
    if missing:
        raise FileNotFoundError('Saved report figures are missing; workbook left unchanged: ' + ', '.join(missing))
    files = report.parent / 'files'
    files.mkdir(exist_ok=True)
    backup_dir = files / 'report_backups'
    backup_dir.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='report-repair-', dir=files) as temporary:
        staged = Path(temporary) / report.name
        shutil.copy2(report, staged)
        for formatter, args in jobs:
            formatter(staged, *args)
        backup = backup_dir / f'{report.stem}-{datetime.now():%Y%m%d-%H%M%S-%f}.xlsx'
        shutil.copy2(report, backup)
        os.replace(staged, report)
    print(f'Report figures restored: {report}\nOriginal workbook backup: {backup}')
    return backup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report', type=Path)
    args = parser.parse_args()
    # Check/repair in a separate process before importing openpyxl here.
    subprocess.run([sys.executable, str(Path(__file__).with_name('reportDependencies.py')), '--ensure'], check=True)
    repair_report(args.report)


if __name__ == '__main__':
    main()
