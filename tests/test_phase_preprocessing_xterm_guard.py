import re
from pathlib import Path


def test_xterm_guard_present(phase_root: Path):
    """xterm is Linux-only. Every `xterm` invocation must sit inside an
    `isunix && ~ismac` branch; Windows must use `start "" /MIN cmd /c`
    without referencing xterm or undefined vars.
    """
    import zipfile
    with zipfile.ZipFile(phase_root / "PHASE_Preprocessing.mlapp") as z:
        xml = z.read("matlab/document.xml").decode("utf-8")

    # No dangling references to variables that were never defined.
    assert "path_2_download_cmd" not in xml, \
        "undefined var path_2_download_cmd still present — ispc branch broken"

    # Every xterm mention must be followed (within a reasonable window) by
    # the unix arm, never be the direct body of an `if ispc` with no
    # elseif ispc/elseif isunix/elseif ismac discipline.
    # Stronger check: structurally, the ispc branch for download should
    # invoke `start "" /MIN cmd /c` (Windows background), and the
    # elseif isunix arm should do xterm.
    download_block = re.search(
        r"STEP 1: Images download started.*?if ispc\s+(.+?)\s+elseif isunix",
        xml, re.DOTALL
    )
    assert download_block, "could not locate STEP 1 ispc branch"
    ispc_body = download_block.group(1)
    assert 'start "" /MIN cmd /c' in ispc_body, \
        f"Windows download must use `start /MIN cmd /c`, got:\n{ispc_body}"
    assert "xterm" not in ispc_body, \
        f"xterm is Linux-only; leaked into ispc branch:\n{ispc_body}"
    # Defense in depth: no nested `if isunix` inside `if ispc` (dead code).
    assert "if isunix" not in ispc_body, \
        f"dead-code `if isunix` inside `if ispc`:\n{ispc_body}"
