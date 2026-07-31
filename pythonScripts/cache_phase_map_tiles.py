#!/usr/bin/env python3
"""Cache the Esri basemap tiles requested by the local MATLAB uihtml map."""

import argparse
import json
import os
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.request import Request, urlopen


SERVICES = {
    "imagery": "https://services.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile",
    "labels": "https://services.arcgisonline.com/ArcGIS/rest/services/Reference/World_Boundaries_and_Places/MapServer/tile",
}
EXTENSIONS = {"imagery": ".jpg", "labels": ".png"}
MAX_REQUESTS = 120


def normalize_request(raw):
    layer = str(raw.get("layer", ""))
    if layer not in SERVICES:
        raise ValueError(f"Unsupported map layer: {layer}")
    zoom = int(raw.get("z"))
    x = int(raw.get("x"))
    y = int(raw.get("y"))
    if not 0 <= zoom <= 18:
        raise ValueError(f"Invalid map zoom: {zoom}")
    limit = 2**zoom
    if not 0 <= x < limit or not 0 <= y < limit:
        raise ValueError(f"Invalid tile coordinates at zoom {zoom}: {x}, {y}")
    key = f"{layer}/{zoom}/{x}/{y}"
    return {"layer": layer, "z": zoom, "x": x, "y": y, "key": key}


def tile_path(cache_root, tile):
    return (
        Path(cache_root)
        / tile["layer"]
        / str(tile["z"])
        / str(tile["x"])
        / f'{tile["y"]}{EXTENSIONS[tile["layer"]]}'
    )


def valid_image(payload):
    return payload.startswith(b"\xff\xd8\xff") or payload.startswith(b"\x89PNG\r\n\x1a\n")


def cache_one(cache_root, tile, timeout=25):
    destination = tile_path(cache_root, tile)
    if destination.is_file() and destination.stat().st_size > 0:
        return {"key": tile["key"], "status": "cached", "path": str(destination)}

    destination.parent.mkdir(parents=True, exist_ok=True)
    url = f'{SERVICES[tile["layer"]]}/{tile["z"]}/{tile["y"]}/{tile["x"]}'
    request = Request(url, headers={"User-Agent": "PHASE-InSAR/1.0 map tile cache"})
    with urlopen(request, timeout=timeout) as response:
        payload = response.read()
    if not valid_image(payload):
        raise RuntimeError("Map service returned data that is not a PNG or JPEG image.")

    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=destination.parent, prefix=f".{destination.name}.", delete=False
        ) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
        os.replace(temporary, destination)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return {"key": tile["key"], "status": "downloaded", "path": str(destination)}


def cache_tiles(cache_root, raw_requests, workers=8):
    if len(raw_requests) > MAX_REQUESTS:
        raise ValueError(f"At most {MAX_REQUESTS} map tiles can be requested at once.")
    unique = {}
    for raw in raw_requests:
        tile = normalize_request(raw)
        unique[tile["key"]] = tile

    successful = []
    failed = []
    with ThreadPoolExecutor(max_workers=max(1, min(int(workers), 12))) as pool:
        futures = {
            pool.submit(cache_one, cache_root, tile): tile for tile in unique.values()
        }
        for future in as_completed(futures):
            tile = futures[future]
            try:
                successful.append(future.result())
            except Exception as exc:
                failed.append({"key": tile["key"], "message": str(exc)})
    successful.sort(key=lambda item: item["key"])
    failed.sort(key=lambda item: item["key"])
    return {"successful": successful, "failed": failed}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-root", required=True)
    parser.add_argument("--requests", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args(argv)

    with open(args.requests, "r", encoding="utf-8") as stream:
        data = json.load(stream)
    raw_requests = data.get("requests", data) if isinstance(data, dict) else data
    if not isinstance(raw_requests, list):
        raise ValueError("Tile request input must be a JSON list.")
    result = cache_tiles(args.cache_root, raw_requests, args.workers)
    with open(args.output, "w", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)


if __name__ == "__main__":
    main()
