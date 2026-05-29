"""Tests for tools/audit_train_windows_compat.py."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

# Make the tool importable. tools/ is not a package.
TOOLS_DIR = Path(__file__).resolve().parent.parent / "tools"
sys.path.insert(0, str(TOOLS_DIR))

import audit_train_windows_compat as audit


def test_strip_removes_line_comment():
    src = "x = 1; % system('rm -rf /')\nfoo();"
    out = audit.strip_matlab_noise(src)
    # Line offsets preserved (newline kept, comment chars become spaces).
    assert "system" not in out
    assert "rm -rf" not in out
    assert out.count("\n") == src.count("\n")


def test_strip_removes_single_quoted_string():
    src = "msg = 'system(rm /tmp/x)'; bar();"
    out = audit.strip_matlab_noise(src)
    # The literal content of the string must not survive.
    assert "system(rm" not in out
    assert "bar" in out  # code outside the string is untouched.


def test_strip_preserves_line_columns():
    src = "abc % comment\nxyz"
    out = audit.strip_matlab_noise(src)
    # Length per line preserved so line numbers and column offsets match.
    assert len(out.splitlines()[0]) == len("abc % comment")
    assert out.splitlines()[1] == "xyz"


def test_strip_handles_escaped_quote_in_string():
    # MATLAB escapes a single quote by doubling it: 'it''s'
    src = "s = 'it''s a test'; q();"
    out = audit.strip_matlab_noise(src)
    assert "test" not in out
    assert "q()" in out


def test_strip_keeps_transpose_operator():
    # `A'` is the transpose operator, not a string. Heuristic: a single
    # quote following an alphanumeric/closing-bracket character is transpose.
    src = "y = A' * b;"
    out = audit.strip_matlab_noise(src)
    assert "A'" in out
    assert "b" in out


def test_scan_detects_system_call(tmp_path: Path):
    f = tmp_path / "x.m"
    f.write_text("x = 1;\nsystem('echo hi');\nz = 2;\n")
    findings = audit.scan_linux_patterns(f)
    assert len(findings) == 1
    assert findings[0].line == 2
    assert findings[0].severity == "HIGH"
    assert findings[0].pattern == "system_call"


def test_scan_detects_aps_systemcall_wrapper(tmp_path: Path):
    f = tmp_path / "x.m"
    f.write_text("aps_systemcall(cmd);\n")
    findings = audit.scan_linux_patterns(f)
    assert len(findings) == 1
    assert findings[0].severity == "HIGH"
    assert findings[0].pattern == "aps_systemcall"


def test_scan_ignores_pattern_in_comment(tmp_path: Path):
    f = tmp_path / "x.m"
    f.write_text("% system('rm /')\nx = 1;\n")
    findings = audit.scan_linux_patterns(f)
    assert findings == []


def test_scan_ignores_pattern_in_string(tmp_path: Path):
    f = tmp_path / "x.m"
    f.write_text("error('use system() to fix');\n")
    # The 'system()' inside the quoted string is stripped.
    # The standalone 'error(' is a builtin, not in our pattern list.
    findings = audit.scan_linux_patterns(f)
    assert findings == []


def test_scan_detects_linux_path_inside_string(tmp_path: Path):
    f = tmp_path / "x.m"
    f.write_text("dest = '/tmp/foo';\n")
    # MATLAB paths ALWAYS live inside string literals — the audit must
    # flag them on the raw source, not the stripped one.
    findings = audit.scan_linux_patterns(f)
    pats = {fi.pattern for fi in findings}
    assert "linux_path_tmp" in pats
    assert all(fi.severity == "MEDIUM" for fi in findings
               if fi.pattern == "linux_path_tmp")


def test_scan_detects_home_expansion_in_fopen(tmp_path: Path):
    # Canonical Windows-blocker pattern from real TRAIN: aps_merra_files.m:50
    # uses `fopen('~/.merrapass','r')`. Must be flagged as MEDIUM.
    f = tmp_path / "x.m"
    f.write_text("fid = fopen('~/.merrapass', 'r');\n")
    findings = audit.scan_linux_patterns(f)
    pats = {fi.pattern for fi in findings}
    assert "home_expansion" in pats


def test_scan_does_not_flag_system_call_inside_string(tmp_path: Path):
    # Regression guard: shell escapes (HIGH) are still stripped-only — they
    # must NOT fire on text that lives inside a string literal.
    f = tmp_path / "x.m"
    f.write_text("error('use system() to fix');\n")
    findings = audit.scan_linux_patterns(f)
    assert all(fi.pattern != "system_call" for fi in findings)


def test_scan_detects_low_severity_shell_command(tmp_path: Path):
    f = tmp_path / "x.m"
    # bare 'wget' identifier outside any string. In real TRAIN this only
    # appears inside system() args, but we want a standalone catch too.
    f.write_text("cmd = wget;\n")  # no string literal here.
    findings = audit.scan_linux_patterns(f)
    assert any(fi.severity == "LOW" and fi.pattern == "shell_wget"
               for fi in findings)


def test_scan_emits_warning_when_eval_or_feval(tmp_path: Path):
    f = tmp_path / "x.m"
    f.write_text("h = feval(name, x);\n")
    findings = audit.scan_linux_patterns(f)
    # `feval` is not a "Linux pattern" but emits a WARNING-severity finding
    # so the report flags incomplete closure for that branch.
    assert any(fi.severity == "WARNING" and fi.pattern == "dynamic_dispatch"
               for fi in findings)


def test_finding_carries_snippet_with_context(tmp_path: Path):
    f = tmp_path / "x.m"
    f.write_text("a = 1;\nb = 2;\nsystem('x');\nc = 3;\nd = 4;\n")
    findings = audit.scan_linux_patterns(f)
    assert len(findings) == 1
    snippet = findings[0].snippet
    # Snippet contains the matched line plus 1 line above and 1 below.
    assert "b = 2;" in snippet
    assert "system" in snippet
    assert "c = 3;" in snippet


def test_index_train_files(tmp_path: Path):
    matlab = tmp_path / "matlab"
    matlab.mkdir()
    (matlab / "foo.m").write_text("function foo(); end")
    (matlab / "bar.m").write_text("function bar(); end")
    (matlab / "sub").mkdir()
    (matlab / "sub" / "baz.m").write_text("function baz(); end")
    idx = audit.index_train_files(tmp_path)
    assert set(idx.keys()) == {"foo", "bar", "baz"}
    assert idx["foo"] == matlab / "foo.m"
    assert idx["baz"] == matlab / "sub" / "baz.m"


def test_index_missing_matlab_dir(tmp_path: Path):
    # No matlab/ subdir → empty index, no exception.
    idx = audit.index_train_files(tmp_path)
    assert idx == {}


def test_extract_calls_finds_function_invocations():
    src = "x = foo(1, 2);\nbar();\nbaz_thing(y);\n"
    calls = audit.extract_calls(src)
    assert "foo" in calls
    assert "bar" in calls
    assert "baz_thing" in calls


def test_extract_calls_ignores_keywords_and_struct_access():
    # `if (...)`, `for (...)`, `while (...)` should not be treated as calls.
    # `obj.method(...)` — we DO want `method` only if it resolves to an .m
    # later, so extract_calls returns it; the index filter handles selection.
    src = "if (x > 0)\n  y = obj.method(x);\nend\nfor i = 1:10\nend\n"
    calls = audit.extract_calls(src)
    assert "if" not in calls
    assert "for" not in calls
    assert "end" not in calls
    assert "method" in calls


def test_extract_calls_strips_strings_and_comments_first():
    src = "% foo()\nbar();\nx = 'baz()';\n"
    calls = audit.extract_calls(src)
    assert "foo" not in calls
    assert "baz" not in calls
    assert "bar" in calls


@pytest.fixture
def phase_root() -> Path:
    return Path("F:/phase")


def _make_train(tmp_path: Path, files: dict[str, str]) -> Path:
    """Helper: create a fake TRAIN tree with given .m contents."""
    matlab = tmp_path / "matlab"
    matlab.mkdir()
    for name, content in files.items():
        (matlab / f"{name}.m").write_text(content)
    return tmp_path


def test_call_graph_single_entry_no_calls(tmp_path: Path):
    root = _make_train(tmp_path, {"aps_linear": "function aps_linear()\nend\n"})
    graph = audit.build_call_graph(root, ["aps_linear"])
    assert len(graph) == 1
    assert next(iter(graph.values())) == set()


def test_call_graph_follows_callees(tmp_path: Path):
    root = _make_train(tmp_path, {
        "aps_linear": "function aps_linear()\n  helper(x);\nend\n",
        "helper": "function helper(x)\n  inner(x);\nend\n",
        "inner": "function inner(x)\nend\n",
        "unreachable": "function unreachable()\nend\n",
    })
    graph = audit.build_call_graph(root, ["aps_linear"])
    names = {p.stem for p in graph}
    assert names == {"aps_linear", "helper", "inner"}
    assert "unreachable" not in names


def test_call_graph_handles_cycles(tmp_path: Path):
    root = _make_train(tmp_path, {
        "aps_linear": "function aps_linear()\n  a();\nend\n",
        "a": "function a()\n  b();\nend\n",
        "b": "function b()\n  a();\nend\n",
    })
    graph = audit.build_call_graph(root, ["aps_linear"])
    assert {p.stem for p in graph} == {"aps_linear", "a", "b"}


def test_call_graph_ignores_unresolved_calls(tmp_path: Path):
    # `sprintf`, `disp` etc. are MATLAB builtins, not in the index.
    root = _make_train(tmp_path, {
        "aps_linear": "function aps_linear()\n  disp('hi');\n  sprintf('x');\nend\n",
    })
    graph = audit.build_call_graph(root, ["aps_linear"])
    assert {p.stem for p in graph} == {"aps_linear"}
    assert next(iter(graph.values())) == set()


def test_call_graph_missing_entry_point_raises(tmp_path: Path):
    root = _make_train(tmp_path, {"foo": "function foo()\nend\n"})
    with pytest.raises(SystemExit) as exc:
        audit.build_call_graph(root, ["aps_linear"])
    assert "aps_linear" in str(exc.value)


def test_call_graph_real_train(phase_root: Path):
    """End-to-end against the real TRAIN clone at F:/phase/TRAIN."""
    train_root = phase_root / "TRAIN"
    if not (train_root / "matlab" / "aps_linear.m").exists():
        pytest.skip("TRAIN clone not present at F:/phase/TRAIN")
    graph = audit.build_call_graph(
        train_root,
        ["aps_linear", "aps_weather_model", "setparm_aps"],
    )
    names = {p.stem for p in graph}
    # Sanity: entry points all reached.
    assert {"aps_linear", "aps_weather_model", "setparm_aps"} <= names
    # Sanity: aps_systemcall is reached (we know it's called transitively).
    assert "aps_systemcall" in names
    # Sanity: closure is non-trivial but bounded by total .m count (88).
    assert 5 <= len(graph) <= 88


def _git_init_with_commit(repo: Path, content: str = "x") -> str:
    """Init a tiny git repo, commit, return the SHA."""
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "config", "user.email", "t@t"], cwd=repo, check=True)
    subprocess.run(["git", "config", "user.name", "t"], cwd=repo, check=True)
    (repo / "f").write_text(content)
    subprocess.run(["git", "add", "f"], cwd=repo, check=True)
    subprocess.run(["git", "commit", "-q", "-m", "x"], cwd=repo, check=True)
    sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo, check=True, capture_output=True, text=True,
    ).stdout.strip()
    return sha


def test_verify_clone_passes_when_head_matches(tmp_path: Path):
    sha = _git_init_with_commit(tmp_path)
    # No exception on match. Pass the full SHA; the function accepts shorts too.
    audit.verify_train_clone(tmp_path, sha)


def test_verify_clone_passes_with_short_sha(tmp_path: Path):
    sha = _git_init_with_commit(tmp_path)
    audit.verify_train_clone(tmp_path, sha[:7])


def test_verify_clone_raises_on_mismatch(tmp_path: Path):
    _git_init_with_commit(tmp_path)
    with pytest.raises(SystemExit) as exc:
        audit.verify_train_clone(tmp_path, "deadbeef")
    assert "deadbeef" in str(exc.value)


def test_verify_clone_raises_when_not_a_repo(tmp_path: Path):
    # tmp_path exists but is not a git repo.
    with pytest.raises(SystemExit):
        audit.verify_train_clone(tmp_path, "abc123")


def test_clone_train_skips_when_dest_exists(tmp_path: Path):
    sha = _git_init_with_commit(tmp_path)
    calls: list[list[str]] = []
    def fake_runner(cmd: list[str]) -> None:
        calls.append(cmd)
    audit.clone_train(tmp_path, sha, runner=fake_runner)
    # Dest exists → no clone, no checkout. Function should return silently
    # after verifying. (verify_train_clone is called internally.)
    assert calls == []


def test_clone_train_invokes_git_when_dest_missing(tmp_path: Path):
    dest = tmp_path / "fresh"
    calls: list[list[str]] = []
    def fake_runner(cmd: list[str]) -> None:
        calls.append(cmd)
        # Simulate clone creating the dir on first call.
        if "clone" in cmd:
            dest.mkdir()
    audit.clone_train(dest, "6c93feb", runner=fake_runner,
                      verify=lambda _d, _c: None)  # skip verify in unit test.
    # Two calls: clone + checkout.
    assert any("clone" in c for c in calls)
    assert any("checkout" in c for c in calls)


def test_render_report_basic(tmp_path: Path):
    findings = [
        audit.Finding(file=tmp_path / "matlab" / "x.m", line=10,
                      pattern="system_call", severity="HIGH",
                      snippet="a\nsystem('x')\nb"),
        audit.Finding(file=tmp_path / "matlab" / "y.m", line=20,
                      pattern="linux_path_tmp", severity="MEDIUM",
                      snippet="c\nx = '/tmp/...'\nd"),
    ]
    md = audit.render_report(
        findings,
        train_commit="6c93feb",
        files_scanned=42,
        train_root=tmp_path,
    )
    assert "# TRAIN Windows compatibility audit" in md
    assert "6c93feb" in md
    assert "42" in md
    assert "## HIGH (1)" in md
    assert "## MEDIUM (1)" in md
    # File sections.
    assert "x.m" in md
    assert "y.m" in md
    # Snippets fenced.
    assert "```" in md


def test_render_report_no_findings(tmp_path: Path):
    md = audit.render_report(
        [], train_commit="6c93feb", files_scanned=42, train_root=tmp_path,
    )
    assert "no findings" in md.lower() or "0 findings" in md.lower()


def test_render_report_paths_relative_to_train_root(tmp_path: Path):
    findings = [
        audit.Finding(file=tmp_path / "matlab" / "deep" / "x.m", line=1,
                      pattern="system_call", severity="HIGH", snippet="x"),
    ]
    md = audit.render_report(
        findings, train_commit="6c93feb", files_scanned=1, train_root=tmp_path,
    )
    # Path printed relative for readability.
    assert "matlab/deep/x.m" in md or "matlab\\deep\\x.m" in md


def test_render_report_includes_limits_section(tmp_path: Path):
    md = audit.render_report(
        [], train_commit="6c93feb", files_scanned=1, train_root=tmp_path,
    )
    assert "Limits of this analysis" in md
    assert "feval" in md.lower() or "dynamic dispatch" in md.lower()


def test_main_writes_report_and_returns_zero_on_clean_train(
    tmp_path: Path, monkeypatch,
):
    # Build a fake TRAIN with one clean entry-point.
    sha = "abcd1234"
    matlab = tmp_path / "matlab"
    matlab.mkdir()
    (matlab / "aps_linear.m").write_text("function aps_linear(); end\n")
    (matlab / "aps_weather_model.m").write_text("function aps_weather_model(); end\n")
    (matlab / "setparm_aps.m").write_text("function setparm_aps(); end\n")
    # Monkeypatch verify to bypass git checks for the fake repo.
    monkeypatch.setattr(audit, "verify_train_clone", lambda _d, _c: None)
    out = tmp_path / "report.md"
    rc = audit.main([
        "--train-path", str(tmp_path),
        "--commit", sha,
        "--output", str(out),
    ])
    assert rc == 0
    assert out.exists()
    text = out.read_text()
    assert "0 findings" in text or "no findings" in text.lower()


def test_main_returns_one_on_high_finding(tmp_path: Path, monkeypatch):
    matlab = tmp_path / "matlab"
    matlab.mkdir()
    (matlab / "aps_linear.m").write_text("function aps_linear()\nsystem('x');\nend\n")
    (matlab / "aps_weather_model.m").write_text("function aps_weather_model(); end\n")
    (matlab / "setparm_aps.m").write_text("function setparm_aps(); end\n")
    monkeypatch.setattr(audit, "verify_train_clone", lambda _d, _c: None)
    out = tmp_path / "report.md"
    rc = audit.main([
        "--train-path", str(tmp_path),
        "--commit", "abc",
        "--output", str(out),
    ])
    assert rc == 1
    text = out.read_text()
    assert "HIGH (1)" in text
    assert "system_call" in text


def test_main_default_output_path(tmp_path: Path, monkeypatch):
    # When --output is omitted, file is written under <phase_root>/reports/.
    matlab = tmp_path / "matlab"
    matlab.mkdir()
    for n in ("aps_linear", "aps_weather_model", "setparm_aps"):
        (matlab / f"{n}.m").write_text(f"function {n}(); end\n")
    monkeypatch.setattr(audit, "verify_train_clone", lambda _d, _c: None)
    monkeypatch.chdir(tmp_path)  # so reports/ is created under tmp_path.
    rc = audit.main([
        "--train-path", str(tmp_path),
        "--commit", "deadbeef",
    ])
    assert rc == 0
    reports = list((tmp_path / "reports").glob("train_audit_*.md"))
    assert len(reports) == 1
    assert "deadbeef"[:7] in reports[0].name
