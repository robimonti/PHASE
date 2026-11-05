"""Exercise actual CSK selection/subset scripts with a controlled SNAP process."""
from pathlib import Path
import os
import subprocess
import sys
import pytest


WRAPPER = r'''
import io,json,runpy,sys,urllib.request
from pathlib import Path
target=Path(sys.argv[1]); config=sys.argv[2]; stub=Path(sys.argv[3])
sys.path.insert(0,str(target.parent))
urllib.request.urlopen=lambda *a,**k:io.BytesIO(json.dumps({'daily':{'precipitation_sum':[0.0]}}).encode())
import phase_subprocess
original=phase_subprocess._ORIGINAL_POPEN
phase_subprocess._ORIGINAL_POPEN=lambda args,**kw:original([sys.executable,str(stub),*args[1:]],**kw)
sys.argv=[str(target),config]
runpy.run_path(str(target),run_name='__main__')
'''


def make_project(tmp_path, phase_root, auto=False):
    project = tmp_path / 'CSK project & data'
    (project / 'slaves').mkdir(parents=True)
    for date in ['20210801', '20210815', '20210901']:
        (project / 'slaves' / f'CSKS1_SCS_B_{date}123456_{date}123500.h5').write_bytes(b'HDF fixture')
    config = project / 'project_master.conf'
    values = {'PROJECTFOLDER': project, 'GRAPHSFOLDER': phase_root / 'PHASE_Preprocessing/snap2stamps/graphs',
              'AUTO_MASTER': int(auto), 'MASTER_DATE': '20210815', 'GPTBIN_PATH': sys.executable,
              'LONMIN': 9, 'LONMAX': 10, 'LATMIN': 45, 'LATMAX': 46, 'CPU': 1, 'CACHE': '128M'}
    config.write_text('\n'.join(f'{k}={v}' for k, v in values.items()))
    return project, config


def run(phase_root, config, script, **environment):
    target = phase_root / 'PHASE_Preprocessing/snap2stamps/bin' / script
    stub = phase_root / 'tests/fixtures/csk_gpt_stub.py'
    return subprocess.run([sys.executable, '-c', WRAPPER, str(target), str(config), str(stub)],
                          capture_output=True, text=True, env={**os.environ, **environment})


@pytest.mark.parametrize('auto', [False, True])
def test_selected_master_reaches_real_subset_script(tmp_path, phase_root, auto):
    project, config = make_project(tmp_path, phase_root, auto)
    selected = run(phase_root, config, 'CSK_master_selector.py')
    assert selected.returncode == 0, selected.stdout + selected.stderr
    assert len(list((project / 'master/20210815').glob('*.h5'))) == 1
    assert len(list((project / 'slaves').glob('*.h5'))) == 2
    subset = run(phase_root, config, 'CSK_subset_master.py')
    assert subset.returncode == 0, subset.stdout + subset.stderr
    assert (project / 'master/20210815_sub.dim').is_file()
    log = (project / 'logs/CSK_subsetM_proc_stdout.log').read_text()
    assert 'TEST SNAP diagnostic: master reader started' in log
    assert 'subset complete' in log
    assert (project / 'logs/CSK_master_selection.log').is_file()


def test_previous_flat_manual_layout_is_supported(tmp_path, phase_root):
    project, config = make_project(tmp_path, phase_root)
    (project / 'master').mkdir()
    source = next((project / 'slaves').glob('*20210815*.h5'))
    source.rename(project / 'master' / source.name)
    result = run(phase_root, config, 'CSK_subset_master.py')
    assert result.returncode == 0, result.stdout + result.stderr
    assert (project / 'master/20210815_sub.dim').is_file()


@pytest.mark.parametrize('env,message', [
    ({'PHASE_TEST_GPT_EXIT': '7'}, 'exit 7'),
    ({'PHASE_TEST_GPT_NO_OUTPUT': '1'}, 'did not create'),
])
def test_snap_failure_is_nonzero_and_saved_to_log(tmp_path, phase_root, env, message):
    project, config = make_project(tmp_path, phase_root)
    assert run(phase_root, config, 'CSK_master_selector.py').returncode == 0
    result = run(phase_root, config, 'CSK_subset_master.py', **env)
    assert result.returncode != 0
    log = (project / 'logs/CSK_subsetM_proc_stdout.log').read_text()
    assert message in log
    assert 'TEST SNAP diagnostic: master reader started' in log
    assert message in result.stderr


def test_missing_manual_date_fails_before_snap_and_keeps_images(tmp_path, phase_root):
    project, config = make_project(tmp_path, phase_root)
    config.write_text(config.read_text().replace('MASTER_DATE=20210815', 'MASTER_DATE=19990101'))
    result = run(phase_root, config, 'CSK_master_selector.py')
    assert result.returncode != 0
    assert '19990101' in (project / 'logs/CSK_master_selection.log').read_text()
    assert len(list((project / 'slaves').glob('*.h5'))) == 3
