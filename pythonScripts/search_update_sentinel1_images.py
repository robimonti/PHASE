#!/usr/bin/env python3
import argparse
import json
import os
import sys
from datetime import datetime, timezone


SENTINEL_PLATFORMS = [
    "Sentinel-1A",
    "Sentinel-1B",
    "Sentinel-1C",
    "Sentinel-1D",
]


def _load_asf_search():
    try:
        import asf_search as asf
    except ImportError as exc:
        raise RuntimeError(
            "Python package 'asf_search' is not installed in the selected environment. "
            "Install it with 'python -m pip install asf_search' or select an environment "
            "where it is already available."
        ) from exc
    return asf


def _write_json(path, payload):
    with open(path, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2)


def _scene_datetime(value):
    if value is None:
        return None
    if isinstance(value, datetime):
        return value
    text = str(value).replace("Z", "+00:00")
    try:
        return datetime.fromisoformat(text)
    except ValueError:
        return None


def _scene_date_yyyymmdd(properties):
    start_time = _scene_datetime(properties.get("startTime"))
    if start_time is None:
        return ""
    return start_time.strftime("%Y%m%d")


def _scene_start_text(properties):
    start_time = _scene_datetime(properties.get("startTime"))
    if start_time is None:
        return str(properties.get("startTime", ""))
    if start_time.tzinfo is None:
        start_time = start_time.replace(tzinfo=timezone.utc)
    return start_time.strftime("%Y-%m-%d %H:%M:%S UTC")


def _as_text(value):
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        return ", ".join(str(item) for item in value)
    return str(value)


def _scene_download_url(scene_name, properties):
    url = properties.get("url") or properties.get("downloadUrl")
    if url:
        return _as_text(url)

    platform_folder = {
        "S1A": "SA",
        "S1B": "SB",
        "S1C": "SC",
        "S1D": "SD",
    }.get(scene_name[:3].upper())
    if not platform_folder:
        return ""

    return f"https://datapool.asf.alaska.edu/SLC/{platform_folder}/{scene_name}.zip"


def _reference_parameters(reference_zip):
    asf = _load_asf_search()
    granule_id = os.path.splitext(os.path.basename(reference_zip))[0]
    reference_results = asf.granule_search(granule_id)
    if not reference_results:
        raise RuntimeError(f"Could not find metadata for reference granule {granule_id}.")

    properties = reference_results[0].properties
    flight_direction = properties.get("flightDirection")
    path_number = properties.get("pathNumber")
    frame_number = properties.get("frameNumber")

    missing = [
        name
        for name, value in (
            ("flightDirection", flight_direction),
            ("pathNumber", path_number),
            ("frameNumber", frame_number),
        )
        if value in (None, "")
    ]
    if missing:
        raise RuntimeError(
            f"Missing {', '.join(missing)} in ASF metadata for {granule_id}."
        )

    return {
        "referenceScene": granule_id,
        "flightDirection": flight_direction,
        "pathNumber": path_number,
        "frameNumber": frame_number,
    }


def search(args):
    asf = _load_asf_search()
    params = _reference_parameters(args.reference)
    excluded_dates = {item for item in args.exclude_dates.split(",") if item}

    start_date = datetime.strptime(args.start_date, "%Y%m%d")
    end_date = datetime.strptime(args.end_date, "%Y%m%d")
    if end_date < start_date:
        raise RuntimeError("End date is earlier than start date.")

    search_results = asf.search(
        platform=SENTINEL_PLATFORMS,
        processingLevel=asf.PRODUCT_TYPE.SLC,
        beamMode=asf.BEAMMODE.IW,
        flightDirection=params["flightDirection"],
        relativeOrbit=params["pathNumber"],
        frame=params["frameNumber"],
        start=start_date.strftime("%Y-%m-%dT00:00:00Z"),
        end=end_date.strftime("%Y-%m-%dT23:59:59Z"),
    )

    scenes = []
    for scene in search_results:
        properties = scene.properties
        scene_date = _scene_date_yyyymmdd(properties)
        scene_name = _as_text(properties.get("sceneName"))
        if scene_date in excluded_dates:
            continue

        scenes.append(
            {
                "sceneName": scene_name,
                "date": scene_date,
                "startTime": _scene_start_text(properties),
                "platform": _as_text(properties.get("platform")),
                "polarization": _as_text(properties.get("polarization")),
                "url": _scene_download_url(scene_name, properties),
            }
        )

    scenes.sort(key=lambda item: (item["date"], item["sceneName"]))
    _write_json(args.output, {"params": params, "results": scenes})


def build_parser():
    parser = argparse.ArgumentParser(
        description="Search Sentinel-1 update imagery from ASF."
    )
    parser.add_argument("--reference", required=True)
    parser.add_argument("--start-date", required=True)
    parser.add_argument("--end-date", required=True)
    parser.add_argument("--exclude-dates", default="")
    parser.add_argument("--output", required=True)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    search(args)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
