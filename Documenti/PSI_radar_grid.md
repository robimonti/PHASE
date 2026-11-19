# Radar grid & PS

Open **StaMPS PSI → Radar grid & PS** after steps 1–4. The panel reads existing
patch outputs; viewing it does not require rerunning PSI. The former diagnostic-export toggle
has been removed from Export; map export belongs to this results panel.

- Drag and scroll to navigate; **Fit area** returns to the candidate extent.
  Radar overlays support zoom up to level 23. Satellite tiles stop at level 18
  and are magnified beyond that; this adds no imagery resolution.
- Colour candidates by selection stage or amplitude dispersion. Filter points
  by stage and switch radar-grid lines on or off. Move over a point to inspect
  its D_A, stage and zero-based azimuth/range indices.
- Grid curves use the same SNAP geolocation validated against the patch
  candidates as the diagnostic exporter. Cell boundaries lie half a pixel from
  zero-based pixel centres. If SNAP rasters are unavailable, the existing
  interpolation from candidate coordinates is explicitly identified.
- At wide extents, grid lines are spaced by multiple cells; the current spacing
  is shown. Zoom in for individual cells. Curves sample at most 512 positions
  per line. At most 40,000 candidates are drawn per viewport; any subsampling is
  disclosed and does not change the complete source table.
- **Export maps & CSV** writes both amplitude-dispersion and selection-stage
  maps for the current geographic extent as PNG and editable MATLAB FIG, plus
  the complete candidate CSV and metadata, under the dataset's `EXPORT` folder.
  Managed projects also receive the files in their regular project exports.
  Colour and stage filters affect the interactive view; both exported maps
  include all candidates in their displayed extent. **Open exports** opens the
  native export folder.

The grid describes radar sampling, not physical resolution. “Retained after
weeding” denotes patch-local survivors; patch merging and subsequent StaMPS
operations can change the final PS set. Selection and scientific processing
remain unchanged. Viewing is read-only; diagnostic exports are written only on request
from this panel. Existing diagnostic files are preserved.
