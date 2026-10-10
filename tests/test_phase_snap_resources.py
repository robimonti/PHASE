"""Runtime resource guards for SNAP GPT on supported platforms."""

import importlib.util


def _helper(phase_root):
    path = phase_root / "PHASE_Preprocessing" / "snap2stamps" / "bin" / "phase_subprocess.py"
    spec = importlib.util.spec_from_file_location("phase_subprocess_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_macos_gpt_guard_caps_cache_and_threads(phase_root):
    helper = _helper(phase_root)
    command, notice = helper.macos_gpt_command(
        ["/Applications/esa-snap/bin/gpt", "graph.xml", "-c", "26G", "-q", "8"],
        physical_bytes=8 * 1024 ** 3, heap_bytes=5 * 1024 ** 3)
    assert command == ["/Applications/esa-snap/bin/gpt", "graph.xml", "-c", "512M", "-q", "2", "-x"]
    assert "26G → 512M" in notice


def test_macos_gpt_guard_preserves_smaller_user_cache(phase_root):
    helper = _helper(phase_root)
    command, _ = helper.macos_gpt_command(
        ["gpt", "graph.xml", "-c", "256M", "-q", "1"],
        physical_bytes=8 * 1024 ** 3, heap_bytes=5 * 1024 ** 3)
    assert command[command.index("-c") + 1] == "256M"
    assert command[command.index("-q") + 1] == "1"
    assert command[-1] == "-x"


def test_macos_gpt_guard_does_not_rewrite_other_commands(phase_root):
    helper = _helper(phase_root)
    command = ["python3", "script.py", "-c", "26G", "-q", "8"]
    assert helper.macos_gpt_command(command) == (command, "")


def test_gpt_guard_applies_to_windows_executable(phase_root):
    helper = _helper(phase_root)
    command, notice = helper.safe_gpt_command(
        ["gpt.exe", "graph.xml", "-c", "26G", "-q", "8"],
        physical_bytes=16 * 1024 ** 3, heap_bytes=4 * 1024 ** 3)
    assert command == ["gpt.exe", "graph.xml", "-c", "512M", "-q", "4", "-x"]
    assert "SNAP memory guard" in notice


def test_gpt_guard_applies_to_linux_executable(phase_root):
    helper = _helper(phase_root)
    command, _ = helper.safe_gpt_command(
        ["/opt/esa-snap/bin/gpt", "graph.xml", "-c", "26G", "-q", "8"],
        physical_bytes=32 * 1024 ** 3, heap_bytes=8 * 1024 ** 3)
    assert command == ["/opt/esa-snap/bin/gpt", "graph.xml", "-c", "1024M", "-q", "8", "-x"]
