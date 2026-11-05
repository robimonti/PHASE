"""Future reports must embed PNGs with the selected Python interpreter."""
import importlib
import subprocess
import sys
import zipfile
import pytest


def test_all_report_sheets_embed_images_without_changing_values(tmp_path, phase_root):
    openpyxl = pytest.importorskip('openpyxl')
    image = pytest.importorskip('PIL.Image')
    folder = tmp_path / 'report with spaces è'
    figures = folder / 'figures'
    figures.mkdir(parents=True)
    names = ['AOI_context', 'AOI_PS', 'AOI_meanDispl_RAW', 'AOI_AvgVelocity', 'AOI_uncertainty', 'ExtrDispl_STA1', 'ExtrDispl_STA2']
    for name in names:
        image.new('RGB', (30, 20), 'red').save(figures / (name + '.png'))
    book = openpyxl.Workbook()
    book.remove(book.active)
    sheets = ['General', 'Raw displacement', 'Modelled displacement', 'Uncertainty', 'AOI anomalies', 'Interpolation', 'Modelled PS']
    for name in sheets:
        sheet = book.create_sheet(name)
        sheet['A1'] = name
    for row, identifier in [(6, 'STA1'), (33, 'STA2')]:
        book['Interpolation'].cell(row, 1, 'ID')
        book['Interpolation'].cell(row + 1, 1, identifier)
        book['Interpolation'].cell(row + 1, 6, 12.3456789)
    book['Modelled displacement']['E33'] = 9531449.42763
    report = folder / 'test_out.xlsx'
    book.save(report)
    expected = {s.title: {c.coordinate: c.value for row in s for c in row if c.value is not None} for s in book}
    book.close()
    general = ';'.join(f'{path},{position},{width},{height}' for path, position, width, height in [
        (phase_root / 'PHASE_logo.png', 'F1', 366, 180),
        (figures / 'AOI_context.png', 'E8', 900, 450),
        (figures / 'AOI_PS.png', 'E30', 900, 450),
    ])
    arguments = [
        [general, 'General'],
        [str(figures / 'AOI_meanDispl_RAW.png'), 'Raw displacement'],
        [str(figures / 'AOI_AvgVelocity.png'), 'Modelled displacement'],
        [str(figures / 'AOI_uncertainty.png'), 'Uncertainty'],
        ['AOI anomalies'],
        ['Interpolation', ','.join(str(figures / f'ExtrDispl_{name}.png') for name in ['STA1', 'STA2'])],
        ['Modelled PS'],
    ]
    for number, args in enumerate(arguments, 1):
        result = subprocess.run([sys.executable, str(phase_root / f'pythonScripts/reportSheet{number}.py'), str(report), *args], capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
    book = openpyxl.load_workbook(report)
    assert {s.title: {c.coordinate: c.value for row in s for c in row if c.value is not None} for s in book} == expected
    assert [len(s._images) for s in book] == [3, 1, 1, 1, 0, 2, 0]
    book.close()
    with zipfile.ZipFile(report) as archive:
        assert len([name for name in archive.namelist() if name.startswith('xl/media/')]) == 8
        assert archive.testzip() is None


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
