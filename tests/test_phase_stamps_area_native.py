"""Functional mask check against an explicitly supplied native StaMPS runtime.

PHASE_TEST_STAMPS_ROOT=/path/to/StaMPS python -m pytest -q tests/test_phase_stamps_area_native.py
"""
import os
from pathlib import Path
import struct
import subprocess
import pytest


def test_native_candidate_selection_and_calamp_honour_mask(tmp_path):
    root = os.environ.get('PHASE_TEST_STAMPS_ROOT')
    if not root:
        pytest.skip('Set PHASE_TEST_STAMPS_ROOT to test the actual native binaries')
    extension = '.exe' if os.name == 'nt' else ''
    selection = Path(root) / 'bin' / ('selpsc_patch' + extension)
    calamp = Path(root) / 'bin' / ('calamp' + extension)
    assert selection.is_file() and calamp.is_file()
    mask = bytes([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1])
    (tmp_path / 'area.mask').write_bytes(mask)
    # Inputs have no spaces: native StaMPS parameter files tokenise filenames.
    for name in ['a.rslc', 'b.rslc']:
        values = [v for excluded in mask for v in (1000.0 if excluded else 10.0, 0.0)]
        (tmp_path / name).write_bytes(struct.pack('>24f', *values))
    (tmp_path / 'selection.in').write_text('0.4\n4\na.rslc 1\nb.rslc 1\n')
    (tmp_path / 'calamp.in').write_text('a.rslc\nb.rslc\n')
    (tmp_path / 'patch.in').write_text('1\n4\n1\n3\n')
    def run(command):
        result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
    run([str(selection), 'selection.in', 'patch.in', 'all.ij', 'all.da', 'all.flt', 'f', '1'])
    run([str(selection), 'selection.in', 'patch.in', 'area.ij', 'area.da', 'area.flt', 'f', '1', 'area.mask'])
    full = (tmp_path / 'all.ij').read_text().splitlines()
    subset = (tmp_path / 'area.ij').read_text().splitlines()
    assert len(full) == 12
    assert {tuple(map(int, row.split()[1:])) for row in subset} == {(1, 1), (1, 2), (2, 1), (2, 2)}
    run([str(calamp), 'calamp.in', '4', 'calamp.out', 'f', '1', 'area.mask'])
    means = [float(row.split()[1]) for row in (tmp_path / 'calamp.out').read_text().splitlines()]
    assert means == pytest.approx([10.0, 10.0])


def test_mt_prep_snap_accepts_mask_and_extracts_only_selected_candidates(tmp_path):
    root = os.environ.get('PHASE_TEST_STAMPS_ROOT')
    if not root:
        pytest.skip('Set PHASE_TEST_STAMPS_ROOT to test the pinned mt_prep_snap runtime')
    root = Path(root)
    export = tmp_path / 'snap'
    for folder in ['rslc', 'geo', 'diff0']:
        (export / folder).mkdir(parents=True)
    for date in ['20200101', '20200113']:
        (export / 'rslc' / (date + '.rslc')).write_bytes(struct.pack('>24f', *([10.0, 0.0] * 12)))
    (export / 'rslc/20200101.rslc.par').write_text('range_samples: 4\nazimuth_lines: 3\n')
    for suffix, values in [('lon', [9, 10, 11, 12] * 3), ('lat', [45] * 4 + [46] * 4 + [47] * 4), ('dem.rdc', [100] * 12)]:
        (export / 'geo' / ('20200101.' + suffix)).write_bytes(struct.pack('>12f', *values))
    (export / 'diff0/20200101_20200113.diff').write_bytes(struct.pack('>24f', *([1.0, 0.0] * 12)))
    work = tmp_path / 'psi'; work.mkdir()
    mask = work / 'psi_area.mask'; mask.write_bytes(bytes([1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1]))
    # Only skip the MATLAB initial-parameter script. All native extraction
    # binaries and the installed mt_prep_snap parser/dispatch execute normally.
    wrapper = 'from stamps import mt_prep_snap as m; import sys; m.run_batch=lambda *a,**k:None; sys.exit(m.main(sys.argv[1:]))'
    import sys
    result = subprocess.run([sys.executable, '-c', wrapper, '20200101', str(export), '0.4', '1', '1', '50', '50', str(mask)],
                            cwd=work, capture_output=True, text=True,
                            env={**os.environ, 'STAMPS': str(root), 'PYTHONPATH': str(root / 'python')})
    assert result.returncode == 0, result.stdout + result.stderr
    candidates = (work / 'PATCH_1/pscands.1.ij').read_text().splitlines()
    assert {tuple(map(int, row.split()[1:])) for row in candidates} == {(1, 1), (1, 2), (2, 1), (2, 2)}
    assert (work / 'PATCH_1/pscands.1.ll').is_file()
    assert (work / 'PATCH_1/pscands.1.ph').is_file()
