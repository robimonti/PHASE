import importlib.util

import pytest


def _module(phase_root):
    path = phase_root / "pythonScripts" / "cache_phase_map_tiles.py"
    spec = importlib.util.spec_from_file_location("cache_phase_map_tiles", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _Response:
    def __init__(self, payload):
        self.payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False

    def read(self):
        return self.payload


def test_tile_cache_downloads_once_and_reuses_local_file(phase_root, tmp_path, monkeypatch):
    module = _module(phase_root)
    calls = []

    def fake_urlopen(request, timeout):
        calls.append((request.full_url, timeout))
        return _Response(b"\xff\xd8\xff" + b"phase-map")

    monkeypatch.setattr(module, "urlopen", fake_urlopen)
    request = {"layer": "imagery", "z": 3, "x": 4, "y": 2}
    result = module.cache_tiles(tmp_path, [request, request], workers=2)

    assert not result["failed"]
    assert result["successful"][0]["status"] == "downloaded"
    assert len(calls) == 1
    assert (tmp_path / "imagery" / "3" / "4" / "2.jpg").read_bytes().startswith(b"\xff\xd8\xff")

    cached = module.cache_tiles(tmp_path, [request], workers=2)
    assert cached["successful"][0]["status"] == "cached"
    assert len(calls) == 1


@pytest.mark.parametrize(
    "tile_request",
    [
        {"layer": "streets", "z": 1, "x": 0, "y": 0},
        {"layer": "labels", "z": 19, "x": 0, "y": 0},
        {"layer": "imagery", "z": 2, "x": 4, "y": 0},
    ],
)
def test_tile_cache_rejects_invalid_requests(phase_root, tile_request):
    module = _module(phase_root)
    with pytest.raises(ValueError):
        module.normalize_request(tile_request)
