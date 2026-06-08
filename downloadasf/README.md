# Download tab usage

## Python Requirements

The backend requires Python and the following Python packages:

- pip install requests
- pip install asf_search
- pip install shapely

Depending on the Python installation, the commands may instead need to be:

- python -m pip install requests
- python -m pip install asf_search
- python -m pip install shapely

or:

- py -m pip install requests
- py -m pip install asf_search
- py -m pip install shapely

The packages must be installed in the same Python environment that MATLAB uses through pyenv.

To check which Python environment MATLAB uses:

matlab
pyenv

## Filters

Press **Show Filters** to open the parameter panel.

Available filters include:

- Product/file type
- Beam mode
- Sub type
- Polarization
- Orbit direction
- Start date and End date
- Path and frame, start and end
- Sampling rate

### Reset all filters

Press:

**Reset all filters**

Parameters are reset to default values.

### Sampling Rate

The Sampling Rate parameter allows limits of number of search results

Examples:

- X images per week
- X images per month
- X images per year

The backend handles the logic for determining which images satisfy the sampling requirement.

---

## Selecting an AOI (Area of Interest)

Two AOI tools are available.

### Draw Rectangle

Press:

**Draw Rectangle**

Clears current AOI

Then drag on the map to create a rectangular AOI.

While drawing mode is active:

- Map navigation is disabled
- The rectangle must be completed before continuing

---

### Draw Polygon

Press:

**Draw Polygon**

Clears current AOI

Then click once for each polygon corner.

To finish the polygon:

- Click again on the first point/node

### Clear AOI

Press:

**Clear AOI**

To remove current AOI from the map.

---

## Running an ASF Search

After selecting parameters and an AOI, press:

**Search ASF**

The application will:

1. Save the selected parameters and AOI to:
   `download/search_request.json`
2. Run the ASF search through the backend
3. Save the search result metadata to:
   `download/search_summary.json`

---

## Search Result Popup

After a successful search, a popup window appears showing:

- Number of SAR images found
- Estimated total download size

Available options:

- Preview Results
- OK

Preview results opens the preview list with the results.
Pressing OK returns to the main page.

---

# Preview Panel

The preview panel can be opened using:

**Show Preview**

The panel contains:

- A table of all search results
- Image selection checkboxes
- Download buttons
- Recommended path/frame information
- Download progress bar

After a download has started, a loading/progress bar appears above the download buttons.

The progress bar shows:

- Download progress in percent
- Estimated download speed
- Estimated remaining download time

---

## Footprints

Press:

**Show Footprints**

to display image footprints on the map.

Press the button again to hide them.

Footprints are automatically shown after a successful search.

---

## Highlighting Footprints

Selecting an image checkbox in the preview table highlights the corresponding footprint on the map.

The highlighted footprint corresponds to the selected image’s:

- Path
- Frame

Only one footprint per unique path/frame combination is shown.

---

## Recommended Path and Frame

The preview panel displays a recommendation in the format:

Recommended: Path X, Frame Y, Direction Z

This recommendation is generated from backend coverage logic and represents the best coverage for the selected AOI.

It is recommended to use these values when selecting images for download.

---

## Login

Before downloading images, the user must log in to ASF.

The login is used to connect the backend API with the ASF services used for downloading SAR images.

Assumes that the user already has an ASF/Earthdata account.

Use the same:

- Username
- Password

that are used for the ASF website.

If:

- One or both input fields are empty, a small error message is displayed
- The username or password is incorrect, the login fails and an error message is shown

The login information is saved to:

`download/login_request.json`

The backend then validates the credentials and returns:

`download/login_result.json`

This file contains information about whether the login attempt was successful or not.

---

# Download Options

Two download options are available:

- Download Selected
- Download All

Download selected only available if one or more images have been selected.

Before downloading, a confirmation popup displays:

- Number of files
- Total download size

---

## Download Restrictions

Images must have the same:

- Path
- Frame
- Flight direction

If incompatible images are selected, the application displays an error message and prevents the download.

---

## Download Storage

Downloaded files are stored in:

`PHASE_Preprocessing/slaves`

Before starting a new download:

- Old irrelevant files may be removed
- Existing duplicate files are preserved

---

## Load Last Download

Press:

Load Last Download

to restore parameters from the most recent successful download.

The following file is used:

`download/last_download_request.json`

This file stores:

- Search parameters
- AOI coordinates

When loaded:

- Filters are restored
- The previous AOI is redrawn on the map

Parameters can then be modified before starting a new search.

---

# Backup JSON Files

After every successful ASF search, backup copies of generated files are created.

These files are currently used internally by the frontend to restore the original full search result after filtering selected images for download.

## search_summary_all.json

Backup of:

`search_summary.json`

Contains the complete search metadata from the latest ASF search.

---

## download_data_all.json

Backup of:

`download_data.json`

Contains the complete download list and URLs returned from the latest ASF search.

---

# Important JSON Files

| File                         | Purpose                                   |
| ---------------------------- | ----------------------------------------- |
| `search_request.json`        | Current ASF search request                |
| `search_summary.json`        | Current active search result metadata     |
| `search_summary_all.json`    | Backup of latest complete search metadata |
| `download_data.json`         | Active download list used by backend      |
| `download_data_all.json`     | Backup of latest complete download list   |
| `last_download_request.json` | Parameters from last successful download  |
| `login_request.json`         | Username and password used for ASF login  |
| `login_result.json`          | Result of backend login validation        |
