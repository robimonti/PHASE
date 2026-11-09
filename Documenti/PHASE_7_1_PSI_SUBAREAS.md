# Optional PSI sub-areas

In **StaMPS PSI → Project**, the **PSI processing area** switch offers:

- **Use full preprocessing area** (default): retain the existing full-area workflow.
- **Select sub-area on map**: draw a polygon on the satellite map over the preview of SNAP radar coverage. Finish drawing before Save or Start. You can drag vertices to refine the boundary, clear it, or fit the map to the dataset coverage.

For a new area, enable **Prepare data for StaMPS** and start at **Step 1**. Start automatically saves the configuration. The source preprocessing stack is read only: PHASE translates the longitude/latitude polygon into a byte mask on the original radar grid, then passes it as the optional eighth argument to `mt_prep_snap` (after the four patch/overlap arguments). StaMPS uses 0 for included pixels and 1 for excluded pixels. Invalid geolocation pixels are excluded. No raster reprojection or spatial resampling is introduced.

Each polygon/master selection runs in an isolated sibling dataset with an `_AOI_<identifier>` suffix, listed in the hub's dataset selector. Different polygons have independent PATCH directories and exports. Returning to full-area mode selects the original dataset. It does not delete sub-area results. Save also selects/creates the appropriate dataset; processing begins only on Start.

The processing folder contains `psi_area.mask` and `psi_area.json`. The JSON records the selected polygon, input-grid identity and selected-pixel count, and is copied beside the published result files. Resume requires the original mask, unchanged area/grid metadata and existing preparation. When the preprocessing geolocation changes, restart that area from Step 1 with data preparation enabled.

The map outline is a lightweight preview envelope of valid geolocation samples; it can include gaps or invalid radar pixels. The mask is evaluated against every pixel, rather than against this preview. Parts of a drawn polygon outside the valid dataset contribute no pixels; an entirely empty selection is rejected.

A smaller area changes the available PS network and potentially the reference, spatial filtering and unwrapping. Its results are not necessarily identical to clipping a full-area PSI result. Choose an area with sufficient coherent points and a suitable stable reference. Masking reduces candidate/PS processing; initial stack reads still use the original raster extent.

All StaMPS parameters are now always shown. Method-dependent choices (for example PERCENT versus DENSITY) still display only their applicable controls.

## Verification

- MATLAB: `phase_stamps_beta.selfTest`, `phase_stamps_beta.psiAreaSelfTest`, `phase_stamps_beta.projectSelfTest`.
- Native runtime: set `PHASE_TEST_STAMPS_ROOT` to the StaMPS runtime and run `python -m pytest -q tests/test_phase_stamps_area_native.py`. Synthetic tests cover the actual `mt_prep_snap` dispatch and native normalization/extraction. Only its MATLAB initial-parameter-script invocation is bypassed in that preparation test.
- Real-image PSI scientific results and the complete Windows/Linux runs still require end-to-end validation.
