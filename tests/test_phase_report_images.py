"""Report images must survive every sheet formatter and safe repair."""
import importlib
from pathlib import Path
import subprocess
import sys
import zipfile
import pytest


def fixture_report(tmp_path, phase_root, monkeypatch):
    openpyxl = pytest.importorskip('openpyxl')
    image = pytest.importorskip('PIL.Image')
    monkeypatch.syspath_prepend(str(phase_root / 'pythonScripts'))
    folder = tmp_path / 'report with spaces, è'
    figures = folder / 'figures'
    figures.mkdir(parents=True)
    for name in ['AOI_context', 'AOI_PS', 'AOI_meanDispl_RAW', 'AOI_AvgVelocity', 'AOI_uncertainty', 'ExtrDispl_STA,1', 'ExtrDispl_STA2']:
        image.new('RGB', (30, 20), 'red').save(figures / (name + '.png'))
    book = openpyxl.Workbook()
    book.remove(book.active)
    for name in ['General', 'Raw displacement', 'Modelled displacement', 'Uncertainty', 'AOI anomalies', 'Interpolation', 'Modelled PS']:
        sheet = book.create_sheet(name)
        sheet['A1'] = name
    for row, identifier in [(6, 'STA,1'), (33, 'STA2')]:
        book['Interpolation'].cell(row, 1, 'ID')
        book['Interpolation'].cell(row + 1, 1, identifier)
        book['Interpolation'].cell(row + 1, 6, 12.3456789)
    book['Modelled displacement']['E33'] = 9531449.42763
    report = folder / 'test_out.xlsx'
    book.save(report)
    return report


def values(book):
    # Formatting may touch blank cells, so compare populated values only.
    return {s.title: {c.coordinate: c.value for row in s for c in row if c.value is not None} for s in book}


def test_image_embedding_and_repair_preserve_values_and_are_repeatable(tmp_path, phase_root, monkeypatch):
    report = fixture_report(tmp_path, phase_root, monkeypatch)
    import openpyxl
    repair = importlib.import_module('repairExcelReport').repair_report
    original = report.read_bytes()
    book = openpyxl.load_workbook(report)
    original_values = values(book)
    book.close()
    backup = repair(report)
    assert backup.read_bytes() == original
    for attempt in range(2):
        book = openpyxl.load_workbook(report)
        assert values(book) == original_values
        assert {s.title: len(s._images) for s in book} == {
            'General': 3, 'Raw displacement': 1, 'Modelled displacement': 1,
            'Uncertainty': 1, 'AOI anomalies': 0, 'Interpolation': 2, 'Modelled PS': 0,
        }
        book.close()
        with zipfile.ZipFile(report) as archive:
            assert len([n for n in archive.namelist() if n.startswith('xl/media/')]) == 8
            assert archive.testzip() is None
        if attempt == 0:
            repair(report)


def test_missing_saved_figure_leaves_report_untouched(tmp_path, phase_root, monkeypatch):
    report = fixture_report(tmp_path, phase_root, monkeypatch)
    original = report.read_bytes()
    (report.parent / 'figures' / 'AOI_uncertainty.png').unlink()
    repair = importlib.import_module('repairExcelReport').repair_report
    with pytest.raises(FileNotFoundError, match='workbook left unchanged'):
        repair(report)
    assert report.read_bytes() == original


def test_dependency_probe_checks_real_png_embedding(phase_root):
    pytest.importorskip('openpyxl')
    pytest.importorskip('PIL.Image')
    result = subprocess.run([sys.executable, str(phase_root / 'pythonScripts/reportDependencies.py')], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert 'image embedding verified' in result.stdout


def test_missing_pillow_is_identified_before_loading_report(monkeypatch, phase_root):
    monkeypatch.syspath_prepend(str(phase_root / 'pythonScripts'))
    dependencies = importlib.import_module('reportDependencies')
    original = dependencies.importlib.import_module

    def without_pillow(name):
        if name == 'PIL.Image':
            raise ImportError('Pillow not installed')
        return original(name)

    monkeypatch.setattr(dependencies.importlib, 'import_module', without_pillow)
    assert 'Pillow' in dependencies.missing_packages()


def test_formatter_failure_during_repair_keeps_original(tmp_path, phase_root, monkeypatch):
    report = fixture_report(tmp_path, phase_root, monkeypatch)
    original = report.read_bytes()
    sheet3 = importlib.import_module('reportSheet3')

    def fail(*args):
        raise RuntimeError('simulated formatter failure')

    monkeypatch.setattr(sheet3, 'format_excel_sheet3', fail)
    repair = importlib.import_module('repairExcelReport').repair_report
    with pytest.raises(RuntimeError, match='simulated formatter failure'):
        repair(report)
    assert report.read_bytes() == original


def test_repair_command_uses_fresh_report_interpreter(tmp_path, phase_root, monkeypatch):
    report = fixture_report(tmp_path, phase_root, monkeypatch)
    result = subprocess.run([sys.executable, str(phase_root / 'pythonScripts/repairExcelReport.py'), str(report)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'image embedding verified' in result.stdout
    assert list((report.parent / 'files/report_backups').glob('*.xlsx'))
