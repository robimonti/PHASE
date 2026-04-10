### Python script to auto-select the optimal CSK Master image
# Evaluates temporal baseline, seasonality (Aug/Sept), and precipitation (Open-Meteo API)

import os
import sys
import glob
import shutil
import json
import urllib.request
import re
from datetime import datetime
from statistics import median

inputfile = sys.argv[1]
bar_message = '\n#####################################################################\n'

print(bar_message)
print("################### CSK MASTER AUTO-SELECTION #######################")
print("## INPUT PARAMETERS ##")

try:
    with open(inputfile, 'r') as in_file:
        for line in in_file.readlines():
            if "PROJECTFOLDER" in line:
                PROJECT = line.split('=')[1].strip()
            if "LONMIN" in line:
                LONMIN = float(line.split('=')[1].strip())
            if "LATMIN" in line:
                LATMIN = float(line.split('=')[1].strip())
            if "LONMAX" in line:
                LONMAX = float(line.split('=')[1].strip())
            if "LATMAX" in line:
                LATMAX = float(line.split('=')[1].strip())
            if "AUTO_MASTER" in line:
                AUTO_MASTER = int(line.split('=')[1].strip())
            if "MASTER_DATE" in line:
                MANUAL_MASTER = line.split('=')[1].strip()
except Exception as e:
    print(f"Error reading config: {e}")
    sys.exit(1)

slaves_dir = os.path.join(PROJECT, 'slaves')
master_dir = os.path.join(PROJECT, 'master')

if not os.path.exists(master_dir):
    os.makedirs(master_dir)

# 1. Manual Override Check
if getattr(sys.modules[__name__], 'AUTO_MASTER', 1) == 0:
    print(f"\n- AUTO-SELECT DISABLED. User forced Master Date: {MANUAL_MASTER}")
    target_date = MANUAL_MASTER
    # Recursive search for manual date
    files_to_move = glob.glob(os.path.join(slaves_dir, "**", f"*{target_date}*.h5"), recursive=True)
    for f in files_to_move:
        shutil.move(f, os.path.join(master_dir, os.path.basename(f)))
    print(f"Moved {len(files_to_move)} file(s) to master folder.")
    sys.exit(0)

print("\n- AUTO-SELECT ENABLED. Analyzing CSK stack...")

# 2. Group downloaded acquisitions by date (Recursive search)
h5_files = glob.glob(os.path.join(slaves_dir, "**", "*.h5"), recursive=True)
if not h5_files:
    print("ERROR: No .h5 files found in the /slaves folder. Download may have failed.")
    sys.exit(1)

date_to_files = {}
for f in h5_files:
    filename = os.path.basename(f)
    match = re.search(r'_(\d{8})\d{6}_', filename)
    if match:
        date_str = match.group(1)
        if date_str not in date_to_files:
            date_to_files[date_str] = []
        date_to_files[date_str].append(f)
    else:
        print(f"Warning: Could not parse date from {filename}")

if not date_to_files:
    sys.exit(1)

date_objects = [datetime.strptime(d, "%Y%m%d") for d in date_to_files.keys()]
date_objects.sort()

# 3. Calculate median date
timestamps = [d.timestamp() for d in date_objects]
median_timestamp = median(timestamps)
median_date = datetime.fromtimestamp(median_timestamp)
print(f"  -> Stack ranges from {date_objects[0].strftime('%Y-%m-%d')} to {date_objects[-1].strftime('%Y-%m-%d')}")
print(f"  -> Ideal temporal center: {median_date.strftime('%Y-%m-%d')}")

# 4. Filter for Late Summer and rank
candidates = []
for d in date_objects:
    days_from_median = abs((d - median_date).days)
    is_late_summer = d.month in [8, 9]
    candidates.append({
        'date_obj': d,
        'date_str': d.strftime("%Y%m%d"),
        'api_str': d.strftime("%Y-%m-%d"),
        'days_from_median': days_from_median,
        'is_late_summer': is_late_summer
    })

candidates.sort(key=lambda x: (not x['is_late_summer'], x['days_from_median']))
top_candidates = candidates[:5]
print(f"\n- Querying Open-Meteo Historical API for top {len(top_candidates)} candidates...")

center_lat = (LATMIN + LATMAX) / 2.0
center_lon = (LONMIN + LONMAX) / 2.0

# 5. Ping API and select optimal date
best_date_str = None
lowest_rain = 999.0

for cand in top_candidates:
    api_url = f"https://archive-api.open-meteo.com/v1/archive?latitude={center_lat}&longitude={center_lon}&start_date={cand['api_str']}&end_date={cand['api_str']}&daily=precipitation_sum&timezone=GMT"

    try:
        req = urllib.request.Request(api_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode())
        rain_mm = data['daily']['precipitation_sum'][0]
        if rain_mm is None: rain_mm = 0.0

        print(f"  -> {cand['api_str']} (Offset: {cand['days_from_median']} days): {rain_mm} mm precipitation.")

        if rain_mm == 0.0:
            best_date_str = cand['date_str']
            lowest_rain = 0.0  # THE FIX: Prevents ghost fallback message
            print(f"\n*** OPTIMAL CSK MASTER FOUND: {best_date_str} (0.0mm rain) ***")
            break
        elif rain_mm < lowest_rain:
            lowest_rain = rain_mm
            best_date_str = cand['date_str']

    except Exception as e:
        print(f"  -> {cand['api_str']}: API Request failed ({e}).")
        if best_date_str is None: best_date_str = cand['date_str']

if best_date_str and best_date_str != candidates[0]['date_str'] and lowest_rain > 0:
    print(f"\n*** SUB-OPTIMAL WEATHER. Selecting driest candidate: {best_date_str} ({lowest_rain}mm rain) ***")

# 6. Move the winning files into a date subfolder in /master
print(f"\n- Moving Master acquisition {best_date_str} to /master folder...")
files_to_move = date_to_files[best_date_str]

# Creates the YYYYMMDD folder inside /master
master_date_folder = os.path.join(master_dir, best_date_str)
if not os.path.exists(master_date_folder):
    os.makedirs(master_date_folder)

for f in files_to_move:
    shutil.move(f, os.path.join(master_date_folder, os.path.basename(f)))

print(f"Successfully staged {len(files_to_move)} file(s) for Master processing.")
print(bar_message)
