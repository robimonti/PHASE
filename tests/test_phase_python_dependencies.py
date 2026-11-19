"""Quality-map support must work for both old environments and new installs."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'PHASE_Preprocessing/snap2stamps/bin/phase_python_dependencies.py'
spec = importlib.util.spec_from_file_location('quality_dependencies', SCRIPT)
deps = importlib.util.module_from_spec(spec)
spec.loader.exec_module(deps)


def test_real_geotiff_writer_and_reader():
    pytest.importorskip('rasterio')
    result = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'write and read verified' in result.stdout


def test_wheel_ignores_foreign_proj_database_in_this_process_only(tmp_path):
    rasterio = pytest.importorskip('rasterio')
    if not (Path(rasterio.__file__).parent / 'proj_data/proj.db').is_file():
        pytest.skip('This system Rasterio build does not bundle PROJ data.')
    env = dict(os.environ, PROJ_DATA=str(tmp_path), PROJ_LIB=str(tmp_path), GDAL_DATA=str(tmp_path))
    result = subprocess.run([sys.executable, str(SCRIPT)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


def test_missing_rasterio_installs_only_in_the_running_interpreter(monkeypatch):
    def missing(name):
        raise ModuleNotFoundError('No module named rasterio', name='rasterio')
    monkeypatch.setattr(deps.importlib, 'import_module', missing)
    calls = []
    monkeypatch.setattr(deps.subprocess, 'run', lambda args, **kw: calls.append(args))
    deps.ensure_geotiff()
    assert calls[0] == [sys.executable, '-m', 'pip', 'install',
                        '--disable-pip-version-check', '--no-input',
                        '--only-binary=:all:', deps.RASTERIO_PACKAGE]
    assert calls[1] == [sys.executable, str(SCRIPT)]


def test_install_failure_stops_before_stack_analysis_with_actionable_error(monkeypatch):
    def missing(name):
        raise ModuleNotFoundError('No module named rasterio', name='rasterio')
    def fail(*args, **kwargs):
        raise subprocess.CalledProcessError(1, args[0])
    monkeypatch.setattr(deps.importlib, 'import_module', missing)
    monkeypatch.setattr(deps.subprocess, 'run', fail)
    with pytest.raises(RuntimeError, match='then retry Build quality map'):
        deps.ensure_geotiff()


def test_installed_binary_failure_is_reported_without_reinstalling(monkeypatch):
    monkeypatch.setattr(deps.importlib, 'import_module', lambda name: None)
    def fail():
        raise OSError('GDAL writer is unavailable')
    monkeypatch.setattr(deps, 'check_geotiff', fail)
    with pytest.raises(RuntimeError, match='failed its write/read check'):
        deps.ensure_geotiff()


@pytest.mark.parametrize('code', [0, 1])
def test_update_prepares_dependencies_and_keeps_optional_failure_recoverable(tmp_path, monkeypatch, code):
    import phase_update
    monkeypatch.setattr(phase_update, 'release_info', lambda: dict(version='v7.2.2', url='fixture', digest='fixture'))
    monkeypatch.setattr(phase_update, 'download_verified', lambda *args: None)
    def extract(archive, destination, tag):
        probe = destination / SCRIPT.relative_to(ROOT)
        probe.parent.mkdir(parents=True)
        probe.write_text('fixture')
    monkeypatch.setattr(phase_update, 'extract_engine', extract)
    calls = []
    def run(args, **kw):
        calls.append(args)
        return SimpleNamespace(returncode=code, stdout='', stderr='fixture missing wheel' if code else '')
    monkeypatch.setattr(phase_update.subprocess, 'run', run)
    result = phase_update.prepare(tmp_path, {'version': 'v7.2.1'})
    assert result['prepared'] and (tmp_path / 'pending-update/engine').is_dir()
    assert calls[0][0] == sys.executable and calls[0][-1] == '--ensure'
    assert bool(result.get('dependencyWarning')) == bool(code)
