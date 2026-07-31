import importlib.util
import json
from types import SimpleNamespace


def _module(phase_root):
    path = phase_root / "pythonScripts" / "phase_download_manager.py"
    spec = importlib.util.spec_from_file_location("phase_download_manager", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _Response:
    def __init__(self, chunks, status=200, headers=None, after_first=None):
        self._chunks = chunks
        self.status_code = status
        self.headers = headers or {}
        self.after_first = after_first

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def raise_for_status(self):
        return None

    def iter_content(self, chunk_size):
        assert chunk_size > 0
        for index, chunk in enumerate(self._chunks):
            yield chunk
            if index == 0 and self.after_first:
                self.after_first()


class _Session:
    def __init__(self, responses):
        self.responses = list(responses)
        self.requests = []

    def get(self, url, **kwargs):
        self.requests.append((url, kwargs))
        return self.responses.pop(0)


def _args(tmp_path, manifest, kind="initial"):
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps({"files": manifest}), encoding="utf-8")
    return SimpleNamespace(
        manifest=str(manifest_path),
        destination=str(tmp_path / "slaves"),
        credentials=str(tmp_path / "credentials.json"),
        progress=str(tmp_path / "progress.json"),
        output=str(tmp_path / "result.json"),
        stop=str(tmp_path / "stop.request"),
        kind=kind,
    )


def test_download_manager_reports_progress_and_preserves_existing_files(phase_root, tmp_path):
    module = _module(phase_root)
    entries = [
        {"name": "existing.zip", "url": "https://example.test/existing.zip", "sizeBytes": 3},
        {"name": "new.zip", "url": "https://example.test/new.zip", "sizeBytes": 6},
    ]
    args = _args(tmp_path, entries)
    destination = tmp_path / "slaves"
    destination.mkdir()
    (destination / "existing.zip").write_bytes(b"old")
    unrelated = destination / "unrelated.zip"
    unrelated.write_bytes(b"keep")
    session = _Session([_Response([b"abc", b"def"], headers={"Content-Length": "6"})])

    result = module.run(args, session=session)

    assert result["status"] == "completed"
    assert result["skipped"] == ["existing.zip"]
    assert result["downloaded"] == ["new.zip"]
    assert unrelated.read_bytes() == b"keep"
    assert (destination / "new.zip").read_bytes() == b"abcdef"
    progress = json.loads((tmp_path / "progress.json").read_text(encoding="utf-8"))
    assert progress["percentage"] == 100
    assert progress["completedFiles"] == 2
    assert progress["totalFiles"] == 2


def test_download_manager_resumes_part_file_with_http_range(phase_root, tmp_path):
    module = _module(phase_root)
    args = _args(
        tmp_path,
        [{"name": "resume.zip", "url": "https://example.test/resume.zip", "sizeBytes": 6}],
        kind="update",
    )
    destination = tmp_path / "slaves"
    destination.mkdir()
    (destination / "resume.zip.part").write_bytes(b"abc")
    session = _Session(
        [_Response([b"def"], status=206, headers={"Content-Range": "bytes 3-5/6"})]
    )

    result = module.run(args, session=session)

    assert result["status"] == "completed"
    assert (destination / "resume.zip").read_bytes() == b"abcdef"
    assert session.requests[0][1]["headers"] == {"Range": "bytes=3-"}


def test_download_manager_honours_stop_request_before_next_file(phase_root, tmp_path):
    module = _module(phase_root)
    args = _args(
        tmp_path,
        [{"name": "stop.zip", "url": "https://example.test/stop.zip", "sizeBytes": 3}],
    )
    (tmp_path / "stop.request").write_text("stop", encoding="utf-8")

    result = module.run(args, session=_Session([]))

    assert result["status"] == "stopped"
    assert result["completedFiles"] == 0
    progress = json.loads((tmp_path / "progress.json").read_text(encoding="utf-8"))
    assert progress["phase"] == "stopped"


def test_download_manager_updates_progress_and_keeps_partial_on_mid_file_stop(
    phase_root, tmp_path
):
    module = _module(phase_root)
    module.PROGRESS_INTERVAL = 0
    args = _args(
        tmp_path,
        [{"name": "large.zip", "url": "https://example.test/large.zip", "sizeBytes": 6}],
    )

    def request_stop():
        progress = json.loads((tmp_path / "progress.json").read_text(encoding="utf-8"))
        assert progress["currentBytes"] == 3
        assert progress["percentage"] == 50
        (tmp_path / "stop.request").write_text("stop", encoding="utf-8")

    response = _Response(
        [b"abc", b"def"],
        headers={"Content-Length": "6"},
        after_first=request_stop,
    )
    result = module.run(args, session=_Session([response]))

    assert result["status"] == "stopped"
    assert (tmp_path / "slaves" / "large.zip.part").read_bytes() == b"abc"
    assert not (tmp_path / "slaves" / "large.zip").exists()
