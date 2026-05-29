"""Static audit of TRAIN's Windows compatibility, scoped to the closure
of TRAIN functions actually invoked by PHASE.

Entry points (verified by extracting PHASE_StaMPS.mlapp):
- aps_linear
- aps_weather_model
- setparm_aps

Run: python tools/audit_train_windows_compat.py [--train-path F:/phase/TRAIN]
"""
from __future__ import annotations

import sys


def strip_matlab_noise(src: str) -> str:
    """Replace MATLAB line comments and single-quoted string literals with
    spaces, preserving line/column offsets.

    Distinguishes the transpose operator (`A'`) from a string opening quote:
    a `'` immediately following an identifier, digit, `)`, `]`, `}`, or `.`
    is treated as transpose. Otherwise it opens a string. MATLAB escapes a
    single quote inside a string by doubling it (`''`).
    """
    out = []
    i = 0
    n = len(src)
    transpose_prev = False  # True if previous non-space char allows transpose.
    while i < n:
        ch = src[i]
        if ch == "%":
            # Line comment to end-of-line.
            while i < n and src[i] != "\n":
                out.append(" ")
                i += 1
            transpose_prev = False
            continue
        if ch == "'" and not transpose_prev:
            # String literal until next un-doubled quote.
            out.append(" ")
            i += 1
            while i < n:
                if src[i] == "'":
                    # Check for doubled quote (escaped).
                    if i + 1 < n and src[i + 1] == "'":
                        out.append("  ")
                        i += 2
                        continue
                    out.append(" ")
                    i += 1
                    break
                if src[i] == "\n":
                    out.append("\n")
                else:
                    out.append(" ")
                i += 1
            transpose_prev = False
            continue
        out.append(ch)
        if ch.isalnum() or ch in ")]}._":
            transpose_prev = True
        elif ch.isspace():
            transpose_prev = False  # space breaks transpose context: `A '` is not transpose.
        else:
            transpose_prev = False
        i += 1
    return "".join(out)


import re
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Finding:
    file: Path
    line: int
    pattern: str
    severity: str  # HIGH | MEDIUM | LOW | WARNING
    snippet: str   # matched line + 1 line of context above and below


# Patterns scanned on the STRIPPED source: shell escapes / MATLAB-builtins
# whose detection would be confused by occurrences inside string literals
# (e.g. error('use system() to fix') must NOT flag system_call).
_PATTERNS_STRIPPED: list[tuple[str, str, re.Pattern[str]]] = [
    # HIGH — direct shell escapes
    ("system_call",       "HIGH",   re.compile(r"\bsystem\s*\(")),
    ("unix_call",         "HIGH",   re.compile(r"\bunix\s*\(")),
    ("popen_call",        "HIGH",   re.compile(r"\bpopen\s*\(")),
    ("aps_systemcall",    "HIGH",   re.compile(r"\baps_systemcall\s*\(")),
    ("bang_shell",        "HIGH",   re.compile(r"(?:^|\s)!\w")),
    ("getenv_home",       "MEDIUM", re.compile(r"\bgetenv\s*\(\s*['\"]HOME['\"]")),
    # WARNING — closure incompleteness signals
    ("dynamic_dispatch",  "WARNING", re.compile(r"\b(?:feval|eval)\s*\(")),
]

# Patterns scanned on the RAW source: Linux paths and shell-command tokens
# that, in MATLAB, ALWAYS live inside string literals (paths passed to
# fopen/system/etc.). Stripping strings would silently hide them — see
# `aps_merra_files.m` `fopen('~/.merrapass')` as a canonical example.
_PATTERNS_RAW: list[tuple[str, str, re.Pattern[str]]] = [
    # MEDIUM — Linux-specific assumptions (path strings passed to MATLAB IO)
    ("linux_path_tmp",    "MEDIUM", re.compile(r"/tmp/")),
    ("linux_path_var",    "MEDIUM", re.compile(r"/var/")),
    ("linux_path_home",   "MEDIUM", re.compile(r"/home/")),
    ("linux_path_usr",    "MEDIUM", re.compile(r"/usr/")),
    ("home_expansion",    "MEDIUM", re.compile(r"~/")),
    ("shell_script_ref",  "MEDIUM", re.compile(r"\.sh\b")),
    # LOW — shell command identifiers (typically inside system(' ... ') args)
    ("shell_wget",        "LOW",    re.compile(r"\bwget\b")),
    ("shell_curl",        "LOW",    re.compile(r"\bcurl\b")),
    ("shell_tar",         "LOW",    re.compile(r"\btar\b")),
    ("shell_gunzip",      "LOW",    re.compile(r"\bgunzip\b")),
    ("shell_rm_flag",     "LOW",    re.compile(r"\brm\s+-")),
    ("shell_cp_flag",     "LOW",    re.compile(r"\bcp\s+-")),
    ("shell_mv",          "LOW",    re.compile(r"\bmv\s+")),
    ("shell_sed",         "LOW",    re.compile(r"\bsed\b")),
    ("shell_awk",         "LOW",    re.compile(r"\bawk\b")),
    ("shell_ncdump",      "LOW",    re.compile(r"\bncdump\b")),
    ("shell_grib",        "LOW",    re.compile(r"\bgrib_")),
]


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        print(f"WARN: {path} not decodable as UTF-8, falling back to "
              f"latin-1", file=sys.stderr)
        return path.read_text(encoding="latin-1")


def _snippet(lines: list[str], lineno: int) -> str:
    """1-indexed line + 1 above + 1 below, joined with newlines."""
    start = max(0, lineno - 2)
    end = min(len(lines), lineno + 1)
    return "\n".join(lines[start:end])


def scan_linux_patterns(file_path: Path) -> list[Finding]:
    raw = _read_text(file_path)
    stripped = strip_matlab_noise(raw)
    raw_lines = raw.splitlines()
    stripped_lines = stripped.splitlines()
    findings: list[Finding] = []
    # Stripped pass: shell escapes (false-positive-prone inside strings).
    for lineno, line in enumerate(stripped_lines, start=1):
        for name, sev, rx in _PATTERNS_STRIPPED:
            if rx.search(line):
                findings.append(Finding(
                    file=file_path,
                    line=lineno,
                    pattern=name,
                    severity=sev,
                    snippet=_snippet(raw_lines, lineno),
                ))
    # Raw pass: Linux paths and shell tokens (live inside string literals).
    for lineno, line in enumerate(raw_lines, start=1):
        for name, sev, rx in _PATTERNS_RAW:
            if rx.search(line):
                findings.append(Finding(
                    file=file_path,
                    line=lineno,
                    pattern=name,
                    severity=sev,
                    snippet=_snippet(raw_lines, lineno),
                ))
    return findings


# MATLAB control-flow keywords that look like function calls.
_MATLAB_KEYWORDS: frozenset[str] = frozenset({
    "if", "elseif", "else", "for", "while", "switch", "case", "otherwise",
    "end", "function", "return", "break", "continue", "try", "catch",
    "global", "persistent", "classdef", "properties", "methods", "events",
    "enumeration", "spmd", "parfor",
})

_CALL_RX = re.compile(r"\b([a-zA-Z][a-zA-Z0-9_]*)\s*\(")


def index_train_files(train_root: Path) -> dict[str, Path]:
    """Map basename (no .m) → absolute Path for every .m under train_root/matlab.

    If two `.m` files share a basename (rare but possible — TRAIN sometimes
    has helper scripts in subdirs), the last one wins after sorted iteration
    and a warning is emitted to stderr so the auditor knows a name clash
    occurred. Iteration is sorted so basename resolution is deterministic
    across platforms (NTFS vs ext4 yield different `rglob` orders), which
    keeps audit reports reproducible.
    """
    matlab_dir = train_root / "matlab"
    if not matlab_dir.is_dir():
        return {}
    idx: dict[str, Path] = {}
    for p in sorted(matlab_dir.rglob("*.m")):
        if p.stem in idx:
            print(f"WARN: duplicate basename {p.stem}: "
                  f"{idx[p.stem]} vs {p}", file=sys.stderr)
        idx[p.stem] = p
    return idx


def extract_calls(src: str) -> set[str]:
    """Return the set of identifiers used as `name(` in src, after stripping
    comments/strings and removing MATLAB keywords."""
    cleaned = strip_matlab_noise(src)
    return {
        m.group(1) for m in _CALL_RX.finditer(cleaned)
        if m.group(1) not in _MATLAB_KEYWORDS
    }


import subprocess
from collections import deque


def build_call_graph(
    train_root: Path,
    entry_points: list[str],
) -> dict[Path, set[str]]:
    """BFS from entry_points across the TRAIN .m files. Returns the reachable
    closure as a dict mapping each reached file to the set of resolved callee
    basenames.
    """
    index = index_train_files(train_root)
    missing = [name for name in entry_points if name not in index]
    if missing:
        raise SystemExit(
            f"Entry-point(s) not found in {train_root}/matlab: {missing}. "
            f"Index has {len(index)} files. Closest matches: "
            f"{[n for n in index if any(m[:4] in n for m in missing)][:5]}"
        )

    graph: dict[Path, set[str]] = {}
    visited: set[str] = set()
    queue: deque[str] = deque(entry_points)
    while queue:
        name = queue.popleft()
        if name in visited:
            continue
        visited.add(name)
        path = index[name]
        src = _read_text(path)
        callees = extract_calls(src)
        # Keep only callees that resolve to a TRAIN .m file. Exclude
        # self-edges: the regex matches the function's own declaration line
        # `function ... = name(args)`, which would otherwise show up as a
        # spurious self-call in the recorded graph data. Cycle prevention
        # is already handled by `visited`, so dropping these here only
        # cleans the reported edge set without changing closure membership.
        resolved = {c for c in callees if c in index and index[c] != path}
        graph[path] = resolved
        for c in resolved:
            if c not in visited:
                queue.append(c)
    return graph


def verify_train_clone(dest: Path, expected_commit: str) -> None:
    """Verify that `dest` is a git repo whose HEAD matches `expected_commit`
    (full or short SHA). Raises SystemExit on mismatch or missing repo."""
    try:
        result = subprocess.run(
            ["git", "-C", str(dest), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        raise SystemExit(
            f"TRAIN path {dest} is not a git repo (or git is not on PATH): {e}"
        )
    head = result.stdout.strip()
    if not head.startswith(expected_commit):
        raise SystemExit(
            f"TRAIN HEAD at {dest} is {head[:12]}, expected {expected_commit}. "
            f"Refusing to overwrite. Either checkout the audited commit or "
            f"pass --commit <new-sha>."
        )


from typing import Callable

TRAIN_REPO_URL = "https://github.com/dbekaert/TRAIN.git"


def _default_git_runner(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def clone_train(
    dest: Path,
    commit: str,
    *,
    runner: Callable[[list[str]], None] = _default_git_runner,
    verify: Callable[[Path, str], None] | None = None,
) -> None:
    """Clone TRAIN into `dest` and check out `commit`. If `dest` already
    exists, only verify that its HEAD matches `commit` (no overwrite).
    The `runner` and `verify` parameters exist for tests.

    `verify` defaults to a lazy lookup of `verify_train_clone` so
    monkeypatching the module attribute in tests works as expected.
    """
    if verify is None:
        # Resolve at call time so tests that monkeypatch
        # audit.verify_train_clone are honored.
        verify = verify_train_clone
    if dest.exists():
        verify(dest, commit)
        return
    runner(["git", "clone", TRAIN_REPO_URL, str(dest)])
    runner(["git", "-C", str(dest), "checkout", commit])


from datetime import date

_SEVERITY_ORDER = ("HIGH", "MEDIUM", "LOW", "WARNING")

_LIMITS_SECTION = """\
## Limits of this analysis

This is a static, regex-based audit. The following caveats apply:

- **Dynamic dispatch** (`feval`, `eval`, function handles) is not followed.
  Files reaching such call sites emit a `WARNING` finding so you can
  manually inspect those branches.
- **Path-separator audits** (`\\` vs `/`) are out of scope — too noisy.
- **MATLAB toolbox-specific behavior** is out of scope.
- **String-leaked patterns**: `'...'` literals are stripped, but multi-line
  or concatenated strings may leak into the LOW band. Treat LOW as a hint,
  not a verdict.

For the runtime behavior of any flagged site, install TRAIN and exercise
the code path in MATLAB.
"""


def render_report(
    findings: list[Finding],
    *,
    train_commit: str,
    files_scanned: int,
    train_root: Path,
) -> str:
    today = date.today().isoformat()
    counts = {s: sum(1 for f in findings if f.severity == s)
              for s in _SEVERITY_ORDER}
    head = [
        "# TRAIN Windows compatibility audit",
        "",
        f"- **TRAIN commit:** `{train_commit}`",
        f"- **Audit date:** {today}",
        f"- **Files scanned:** {files_scanned}",
        f"- **Findings:** "
        + ", ".join(f"{s}={counts[s]}" for s in _SEVERITY_ORDER),
        "",
        _LIMITS_SECTION,
    ]
    if not findings:
        head.append("## Result")
        head.append("")
        head.append("0 findings — closure is clean against the audited patterns.")
        return "\n".join(head) + "\n"

    body: list[str] = []
    for sev in _SEVERITY_ORDER:
        bucket = [f for f in findings if f.severity == sev]
        if not bucket:
            continue
        body.append(f"## {sev} ({len(bucket)})")
        body.append("")
        # Group by file for readability.
        by_file: dict[Path, list[Finding]] = {}
        for f in bucket:
            by_file.setdefault(f.file, []).append(f)
        for file_path, items in sorted(by_file.items()):
            try:
                rel = file_path.relative_to(train_root).as_posix()
            except ValueError:
                rel = str(file_path)
            body.append(f"### `{rel}`")
            body.append("")
            for item in items:
                body.append(f"- **line {item.line}** — `{item.pattern}`")
                body.append("")
                body.append("  ```matlab")
                for ln in item.snippet.splitlines():
                    body.append(f"  {ln}")
                body.append("  ```")
                body.append("")
    return "\n".join(head + body) + "\n"


import argparse


ENTRY_POINTS = ["aps_linear", "aps_weather_model", "setparm_aps"]
DEFAULT_TRAIN_PATH = Path("F:/phase/TRAIN")
DEFAULT_COMMIT = "6c93feb"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Audit TRAIN's Windows compatibility, scoped to the "
                    "transitive closure of the TRAIN functions invoked by PHASE."
    )
    parser.add_argument(
        "--train-path", type=Path, default=DEFAULT_TRAIN_PATH,
        help="Path to TRAIN clone (default: %(default)s).",
    )
    parser.add_argument(
        "--commit", default=DEFAULT_COMMIT,
        help="TRAIN commit to audit (default: %(default)s).",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Output markdown path. Default: reports/train_audit_<sha>_<date>.md",
    )
    args = parser.parse_args(argv)

    # Ensure TRAIN is present at the expected commit. Does not clone if dest
    # exists; user must remove it explicitly to re-clone.
    clone_train(args.train_path, args.commit)

    graph = build_call_graph(args.train_path, ENTRY_POINTS)
    findings: list[Finding] = []
    for file_path in graph:
        findings.extend(scan_linux_patterns(file_path))

    md = render_report(
        findings,
        train_commit=args.commit,
        files_scanned=len(graph),
        train_root=args.train_path,
    )

    if args.output is None:
        reports_dir = Path.cwd() / "reports"
        reports_dir.mkdir(exist_ok=True)
        sha7 = args.commit[:7]
        out = reports_dir / f"train_audit_{sha7}_{date.today().isoformat()}.md"
    else:
        out = args.output
        out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(md, encoding="utf-8")

    # Console summary.
    counts = {s: sum(1 for f in findings if f.severity == s)
              for s in _SEVERITY_ORDER}
    print(f"Audit complete. Files scanned: {len(graph)}. "
          f"Findings: " + ", ".join(f"{s}={counts[s]}" for s in _SEVERITY_ORDER))
    # Per-file top offenders (spec: console summary includes per-file counts).
    if findings:
        per_file: dict[Path, int] = {}
        for f in findings:
            per_file[f.file] = per_file.get(f.file, 0) + 1
        top = sorted(per_file.items(), key=lambda kv: -kv[1])[:5]
        print("Top files by finding count:")
        for path, cnt in top:
            try:
                rel = path.relative_to(args.train_path).as_posix()
            except ValueError:
                rel = str(path)
            print(f"  {cnt:>3}  {rel}")
    print(f"Report: {out}")

    return 1 if counts["HIGH"] > 0 else 0


if __name__ == "__main__":
    raise SystemExit(main())
