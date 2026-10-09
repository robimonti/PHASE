# PSI export corrections

PHASE asks for an explicit **standard** or **corrected** displacement export.
The choice affects both the full CSV/XLSX export and the optional TS Points
picker. The actual StaMPS `ps_plot` value type is recorded in
`EXPORT/<name>_series.json`.

| Export choice | Available source | `ps_plot` value type | Export-stage subtraction |
| --- | --- | --- | --- |
| Standard | Step 7 or 8 | `v-do` | DEM/SCLA error and orbital ramp; no extra `a`/`s` subtraction |
| Corrected | TRAIN enabled and available | `v-dao` | TRAIN tropospheric estimate (`a`), plus `d` and `o` |
| Corrected | No TRAIN, Step 8 completed | `v-dso` | Step 8 slave spatially correlated noise (`s`), plus `d` and `o` |

If both TRAIN and Step 8 are configured, corrected export uses `v-dao` only;
the two estimates are not silently stacked. If the requested TRAIN correction
is unavailable, corrected export stops rather than falling back to `v-do`.
Existing configurations without the new choice migrate to their historical
behavior: TRAIN runs use corrected `v-dao`; runs without TRAIN use standard
`v-do`.

The terms are deliberately precise. StaMPS Step 8 filters *spatially
correlated noise*; atmospheric delay is a major component, but the estimate
is not guaranteed to be purely atmospheric and may absorb spatially
correlated deformation. Likewise, standard `v-do` means no additional
atmospheric subtraction by `ps_plot`; it does not guarantee that every
upstream processing operation was atmosphere-free.

The `ph_mm` time series from StaMPS are already in millimetres, and the
`ph_disp` velocity is in millimetres per year. PHASE inserts a zero-valued
master epoch only when it is absent from the StaMPS time axis and checks that
the date, displacement, velocity and coordinate dimensions agree.

Primary references:

- [StaMPS `ps_plot.m`](https://github.com/dbekaert/StaMPS/blob/master/matlab/ps_plot.m), especially the `v-do`, `v-dso`, `v-dao` branches.
- [StaMPS `stamps.m`](https://github.com/dbekaert/StaMPS/blob/master/matlab/stamps.m), which defines Step 8 as spatially correlated noise filtering.
- [StaMPS manual, section 4.8](https://forum.step.esa.int/uploads/default/original/2X/5/5e96ab63d4b9da43fb5de64f3306d8f31606b29e.pdf).
- [StaMPS `ts_flaghelper.m`](https://github.com/dbekaert/StaMPS/blob/master/matlab/ts_flaghelper.m), which converts phase to millimetres for the time-series MAT file.
