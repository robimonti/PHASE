"""Structural checks for TRAIN Windows-port Changes #1–#3 in PHASE_StaMPS.mlapp.

Requires phase_root fixture, defined in tests/conftest.py (existing).
"""
import re
import zipfile
from pathlib import Path


def _read_xml(mlapp: Path) -> str:
    with zipfile.ZipFile(mlapp) as z:
        return z.read("matlab/document.xml").decode("utf-8")


def test_change1_guard_inserted_at_correct_anchor(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    # Guard must sit directly after the load('input_StaMPS.mat', ...) block
    # (i.e. between 'ph_output'); and the 'Conversion of variables' comment).
    anchor = "'ref_radius_w', 'ph_output');"
    guard_start = "% begin TRAIN availability check"
    conversion = "% Conversion of variables' format"
    # Friendly assertions BEFORE index() so pytest reports "Change #1 not
    # applied" rather than a bare ValueError.
    assert anchor in xml, f"load() anchor {anchor!r} not found in mlapp"
    assert guard_start in xml, "Change #1 begin sentinel missing — guard not applied"
    assert conversion in xml, f"Conversion comment {conversion!r} missing"
    # Ordering check: anchor → guard_start → conversion
    i_anchor = xml.index(anchor)
    i_guard = xml.index(guard_start)
    i_conv = xml.index(conversion)
    assert i_anchor < i_guard < i_conv, \
        "Change #1 guard is not between the load() and the Conversion section"


def test_change1_multi_function_probe(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    # All three TRAIN symbols must be probed (defeats false positives)
    assert "~isempty(which('aps_linear'))" in xml
    assert "~isempty(which('aps_weather_model'))" in xml
    assert "~isempty(which('setparm_aps'))" in xml


def test_change1_warning_id_and_sprintf_safe_form(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "'StaMPS:phase:trainNotAvailable'" in xml, "Canonical warning id missing"
    assert "'%s'" in xml, \
        "Change #1 must use warning(id, '%s', msg) to prevent %-mangling"


def test_change1_local_only_degradation(phase_root: Path):
    """Change #1 must NOT touch app state (preserves user intent for re-run).

    Note: scope all assertions to the guard slice. The bare substring
    'train_flag = 1;' already appears elsewhere in the unmodified file
    (e.g. line 1308: `app.train_flag = 1;` in the checkbox callback),
    so a full-file `in xml` check would pass even if Change #1 was
    never applied. We guard against that with index()-based friendly
    assertions before slicing.
    """
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "% begin TRAIN availability check" in xml, \
        "Change #1 not applied: begin sentinel missing"
    assert "% end TRAIN availability check" in xml, \
        "Change #1 not applied: end sentinel missing"
    guard_start = xml.index("% begin TRAIN availability check")
    guard_end = xml.index("% end TRAIN availability check")
    assert guard_end > guard_start, "sentinel ordering inverted"
    guard = xml[guard_start:guard_end]
    # The local downgrade must be inside the guard
    assert "train_flag = 1;" in guard, \
        "Change #1 guard missing the local `train_flag = 1;` downgrade"
    # No app-state mutation (ASSIGNMENT, not just string mention —
    # the guard comments legitimately reference app.train_flag in prose).
    assert not re.search(r"app\.train_flag\s*=", guard), \
        "Change #1 must not assign to app.train_flag (local degradation only)"
    assert not re.search(r"TRAINatmosphericcorrectionCheckBox\.Value\s*=", guard), \
        "Change #1 must not assign to the checkbox value"


def test_change1_sentinel_comment(phase_root: Path):
    """The sentinel is consumed by Tier-2; don't rename without updating Tier-2."""
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "% end TRAIN availability check" in xml


def test_change2_misleading_stub_removed(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "TRAIN is not available on Windows" not in xml
    assert "StaMPS:phase:trainUnavailable" not in xml
    # Structure preservation: the if/else/end block's trailing inline
    # comment must still be there.
    assert "end % source all the softwares and prepare the data" in xml


def test_change3_stale_comment_updated(phase_root: Path):
    xml = _read_xml(phase_root / "PHASE_Preprocessing/PHASE_StaMPS.mlapp")
    assert "TRAIN not ported to Windows" not in xml
    assert "TRAIN on MATLABPATH, no shell config to source" in xml
