# PHASE 7.1 — Custom deterministic models

Automatic is the default. The four original automatic scientific functions retain their native
modelling and outlier-rejection paths; the optional residual check refits the
existing polynomial, including its uncertainty calculation. Each separate `STmodel_*_custom.m` entry copies its original
outlier-rejection block verbatim and then calls the custom two-stage path.

## Workflow

1. PHASE provides PS and query/grid coordinates and acquisition times.
2. After the original outlier rejection, evaluate the assigned deterministic
   model or fit its free parameters once. Parameter optimisation can iterate
   internally; there is no outer model/covariance alternation.
3. Subtract the model. DET applies the original residual spline fit and
   selection; STC applies the original empirical product covariance, analytical
   covariance fits and collocation equations. Restore the deterministic
   contribution at every prediction location.

The native residual implementations are in `nativeDet1D.m`,
`nativeDet2D.m`, `nativeCovariance1D.m` and `nativeCovariance.m`.
Spline selection uses the existing MDL/variance/F/chi-square controls and
regularisation settings. STC uses the existing covariance-bin and
Gaussian/exponential/Gaussian+cosine controls. The earlier experimental
variogram fit and GCV residual spline implementation have been removed.

The native 1D geoSplinter operator is calibrated on the observation-design
column space, with the selected geometry and lambda fixed. This preserves
its actual Laplacian regularisation for covariance propagation. It requires
additional native linear solves (up to the design rank), not a joint
deterministic/covariance optimisation loop. The 2D adapter retains entries
for skipped spline candidates before sorting; automatic source is untouched.
Covariance adapters use serial accumulation, safe plot sampling and positive
distance-based automatic spatial bins; empirical/fitting functions are reused.

## Interface and model contract

In Displacement Modeling choose DET or STC, then
**Model design → Deterministic model → Custom**.
Use an expression and the parameter table, or a self-contained MATLAB file:

```matlab
function displacement = myModel(x,y,s,t,p)
displacement = (p.a+p.b*x).*(-expm1(-max(t-p.te,0)./p.tau));
end
```

Output is one finite real LOS displacement in mm per input row. A selected
MATLAB file executes as user code. Expressions support numerical arithmetic,
comparisons and the listed elementary numerical functions; use element-wise
operators where needed. Do not estimate a time-independent physical offset:
temporal referencing removes it and makes it unidentifiable.

PHASE evaluates the same function at PS, grid and query positions.
For the normal PHASE PSI-to-model workflow, no manual coordinate preparation
is required: use the default metric coordinates, leave the origins empty,
and define the model and its parameters. The optional **Coordinate settings**
section is collapsed by default and opens automatically for non-default values.
Centering is a coordinate convention, not a separate data correction. Free
coefficients are estimated in that convention. Fixed external coefficients
must retain the units and origin for which they were defined; changing an
origin without transforming those coefficients generally changes the model.
Choose project metric coordinates or longitude/latitude in degrees.
Model origins can be assigned; empty origins use the mean selected PS coordinate.
The chainage `s` is available for 1D. Time is days from the first acquisition
(t0); there is no separate time-offset control. PHASE references physical model values
to the input master epoch. Existing workbook output conventions are retained.

**DSC example:** select
`examples/custom_models/DSC_quadratic_example.m`, longitude/latitude,
X origin `152.20890917`, Y origin `-4.23823009`, no parameters (`[]`).
The supplied coefficients used a model time origin of 2021-01-02. For exact
reuse, the input must start on that date; otherwise adapt the polynomial
inside the user function or refit its parameters. PHASE always supplies t=0
at the first input acquisition. Coefficients in this file
are the rounded values supplied in DSC.txt; at the given target and day 540
the absolute model evaluates to approximately 127.962858 mm. PHASE's relative
value also subtracts the model at the input reference epoch. Obtain full
precision coefficients before using the file as an exact reproduction of
the original fit.

Preview is schematic, using current assigned/initial parameter values and
illustrative model coordinates; it does not fit the input data.

## Residual mean

Default **Model only — diagnose residual mean** performs no extra reduction.
Exports give global mean/RMS, raw and centered second moments, mean by
acquisition, and mean by PS. The warning ratio is descriptive, not an
independence-based statistical significance test or a proof about E[r].

Optional **Estimate constant residual offset** adds a nuisance intercept,
jointly estimated with the free deterministic parameters before the residual
stage. It is included in the estimator sensitivity and covariance propagation.
This constant cancels when output displacement differences are formed; its
effect on residual estimation remains accounted for. No epoch-specific
planes or temporal trend are silently fitted.

The covariance adapters retain the original cleanObs policy: spatial
nonzero-lag products use epoch-centered samples, the spatial lag-zero term
uses the uncentered second moment, and temporal products are unchanged.
These choices do not guarantee that an arbitrary custom model has zero-mean
residuals. If the mean varies systematically over space/time, refine the
deterministic model; a constant offset cannot remove that structure.

## Parameter uncertainty

Fixed parameters default to zero assigned covariance. Optional standard
deviations in the table define a diagonal covariance; a supplied symmetric
positive semidefinite JSON matrix overrides it, in table order. Assigned
covariance may involve only fixed parameters and is assumed independent of
input observation errors. This assumption must be checked if the values were
fitted externally from the same observations. PHASE-estimated parameters use
their shared-data estimator sensitivity instead of this independence assumption.

Let L=dp_hat/dy, J=dm_obs/dp, Jq the derivative of the output mean difference
and W the native residual prediction operator. Define F=Jq L,
G=W(I-J L), H=F+G. STC's local error covariance is:

```
Kqq - H Koq - Kqo H' + H Cyy H'
```

DET's sampling covariance is H N H'. The deterministic/residual cross term is
included. For independently assigned parameter covariance P, add
(Jq-W Jobs) P (Jq-W Jobs)' to the total; the decomposition includes the
negative cross term caused by subtracting and restoring those same parameters.
Assigned uncertainties therefore are not simply added twice.
When assigned uncertain and estimated parameters coexist, PHASE also
propagates how the estimated parameters change with the assigned values,
using implicit differentiation of the same observed Hessian. The full
parameter covariance and output covariance include these induced cross terms.

Nonlinear sensitivity uses the observed least-squares Hessian including
residual curvature. Singular, indefinite or boundary solutions are rejected.
Uncertainty is local and conditional on the selected spline geometry,
lambda, covariance models/parameters and outlier mask. It excludes
hyperparameter/model-selection uncertainty, smoothing bias and model
discrepancy. Numerical covariance stabilization is recorded explicitly.
The native covariance/noise assumptions are retained; additional input
reference-area uncertainty is not inferred from the Excel.

## Residual controls (v7.1.3)

The 1D model / 2D model pages expose the same controls for automatic/custom spatio-temporal
DET and STC, in 1D and 2D. Automatic-engine defaults to Auto, with no range constraints.
Custom workflows select their residual field automatically. The original outlier rejection is unchanged.

**Automatic-engine residual temporal trend:** Diagnose only, Auto, Linear or Quadratic.
The check fits the residual displacement series, never the covariance curve.
Auto compares both the AOI mean and individual PS with constant, linear and quadratic descriptions
using contiguous held-out time blocks. A polynomial requires at least 20%
less prediction MSE than a constant; degree two also requires 10% less than
degree one, and a fitted excursion of at least max(0.5 mm, 10% of residual RMS).
These are documented selection heuristics, not iid p-values or a guarantee of
stationarity. Held-out per-PS errors are aggregated before selection so opposite local trends remain visible. Per-PS coefficients and aggregate scores are exported alongside the AOI mean diagnostics.

Custom DET/STC automatically checks an additional residual deterministic field;
there are no user controls for its coefficients or complexity. Candidates use
linear/quadratic time terms with common, affine/quadratic spatial coefficients,
or the existing native cubic spatial spline basis at coarse resolutions. All
terms are referenced to the input master. They are restored at prediction locations by default; the custom panel now allows excluding them from the prediction target.
The fixed custom model is unchanged. Regularisation and complexity are selected
using five contiguous time blocks and five contiguous spatial groups. The held-out temporal error must improve by at least 10% over no additional
correction. The simplest field within 10% of the best temporal score is retained;
spatial holdouts guide regularisation among temporally equivalent fits and
report interpolation/coverage limitations, without vetoing removal of an
observed temporal trend. Tiny effects
and numerically weak fits are rejected. Scores are conditional on the upstream
custom fit and outlier mask, not an independent validation of the full pipeline.

If no spatial field passes temporal validation, the separately supported common temporal mean
is retained; otherwise no additional correction is applied. The reason and all
candidate scores are saved. The method never forces an empirical covariance
tail towards zero by increasing deterministic complexity until it disappears.

For the selected weighted/regularised coefficient map M, the derivative is
`Ltrend = M - M*Jcustom*Lcustom`. The augmented design/derivative enters the
existing prediction equations, retaining cross terms with the custom fit and
collocation/splines, including assigned parameter uncertainty. Uncertainty is
conditional on field selection and regularisation: selection uncertainty and
regularisation/model bias are not included. Observation, custom-model, added
field and residual contributions are exported separately. Mean and covariance
fit checks can request scientific review; they do not certify stationarity.
Extrapolation outside observation coverage still requires care.

Automatic correction reuses the existing polynomial basis and uncertainty
calculation. It optionally refits the selected native model once, retaining
its time or time-squared terms. Local checks can activate time-space interactions up to native total degree three, admitting only independent added columns and retaining the original spatial degree. Corrections are restored by native prediction code. With automatic STC
`residualAtm`, per-epoch atmospheric means and an additional common temporal
trend are not separately identifiable. Auto therefore diagnoses without refitting in that mode; explicit linear/quadratic correction is rejected. Use `cleanObs` when the temporal trend is to be retained as deformation.

**Correlation lengths (STC):** optional suggested spatial/temporal lengths,
refinement percentage (default ±25%) and independent maximum lengths.
Lengths mean the e-folding distance/time of the kernel envelope: L = 1/c2 for
exponential; L = 1/sqrt(c2) for Gaussian and Gaussian-cosine. A maximum is a
hard bound, not a lag-bin size. When enabled, a bounded weighted fit refines
the native kernel against empirical nonzero-lag bins using the native pair
cofactors. The native amplitude/nugget combination is retained. Cosine period
is kept at the native fit. Without these controls the original fit is used
exactly. A constrained range is not evidence of a good fit; inspect the figure
and the saved range/bound diagnostics. Fit/model-selection uncertainty remains
outside the conditional prediction uncertainty.

Background: [NIST stationarity](https://www.itl.nist.gov/div898/handbook/pmc/section4/pmc442.htm)
and [MathWorks mean plus covariance model](https://www.mathworks.com/help/stats/gaussian-process-regression-models.html).

## Additional outputs

Diagnostic data files below are inside `output_NNN/files/customModel/`. All figures
are in the shared `output_NNN/figures/`, as PNG plus editable FIG. The general Excel
report remains at `output_NNN/`, alongside the existing output subfolders.


- `custom_model.json`: model, parameters/covariance, units, origins, references,
  native residual settings and uncertainty scope.
- `custom_model_results.mat`: retained/original observations, components,
  coordinates and variance decomposition.
- `custom_model_components.csv`: component values and signed covariance terms.
- `custom_residual_mean*.csv/json`: residual-mean diagnostics.
- `custom_model_source/`: copied custom MATLAB function, when selected.
- `figures/STC_cov1D.*` or `figures/STC_cov2D.*`: the original native
  two-panel covariance figure (PNG and editable FIG), with empirical, signal
  and noise. `custom_covariance_time.csv` and `custom_covariance_space.csv`
  retain empirical/kernel values and pair counts without altering the figure.
- `custom_observation_components.csv`: observed = custom model + extra temporal
  trend + residual actually used for fitting, for every retained observation.
- `temporal_before*` / `temporal_after*`: trend diagnostics, including per-PS
  coefficients; after files exist only when a correction is applied.
- Residual diagnostic figures are not generated. Numerical mean/trend checks remain in CSV/JSON.
- `custom_query_timeseries.csv`, when query points are enabled: dated point
  predictions with station ID, coordinates, deterministic/residual components
  and conditional standard deviation, referenced to the first input date.

In **Outputs**, enable **Interpolate query-point time series**, then select
**Load TXT / CSV** or **Enter points**. The GUI list accepts ID, longitude and
latitude in WGS84 degrees; Save writes `query_points_gui.txt` in the project
folder. Both sources follow the same prediction and AOI-selection path.
Points outside the AOI follow the existing warning/exclusion rules. For
custom DET/STC these are evaluations at the requested coordinates, not the
nearest output-grid cell. Re-run modelling after changing query points.

For GNSS comparison use the same LOS direction/sign, spatial reference and
temporal zero; the custom query CSV uses the first input acquisition as zero.
An externally fitted deterministic model does not by itself guarantee
agreement after a different observation subset and residual model are used.

Older STC covariance figures can be recovered from their JSON without refitting:

```matlab
phase_model_custom.exportSavedCovariance('custom_model.json', 'diagnostics')
```

Old results without pair counts show that information as unavailable.
Covariance fitting does not enforce zero residual mean; inspect mean-by-date
and mean-by-PS diagnostics separately. Automatic native temporal checks are enabled by default.
Automatic-engine diagnostics are saved in `files/residualDiagnostics/`.

The normal workbook and spatial export pipeline also receive custom results.

## Verification

```matlab
phase_model_science.residualFieldSelfTest
phase_model_science.selfTest
phase_model_science.nativeIntegrationSelfTest
phase_model_custom.scientificSelfTest
phase_model_custom.filterIntegrationSelfTest
phase_model_beta.selfTest
```

Tests cover a covariance matrix oracle and 20,000 Monte Carlo draws, assigned
uncertainty, optional residual mean, bounded nonlinear fitting and the DSC
example. The four native entry integrations exercise outlier filtering plus
native residual modelling. Python tests check automatic source isolation,
verbatim filter copies and reuse of the native residual selection code.

Prediction evaluates the deterministic model, the residual field and its
uncertainty. Native STC uses a dense covariance factorisation, followed by
memory-bounded blocks of space-time queries. Algebraically equivalent
Cholesky-whitened equations avoid repeated full covariance multiplications.
The monitor separates covariance fitting, factorisation, grid predictions,
PS predictions and optional query predictions. Percentages indicate completed
work within each phase, not an estimated fraction of total runtime.

The STC memory budget is specified in GiB. The guard estimates six dense
double-precision observation matrices plus a 128 MiB prediction workspace;
this is a working estimate, not a hard cap on total MATLAB/process memory.
Leave RAM for native covariance fitting, reports and other applications.
AOI size controls retained observation count; grid resolution and the output
time step in days control prediction count. No observations are silently dropped.

The custom path is released for real-project testing. Validate the full
GUI/report path on a real project and review residual-mean diagnostics before
relying on its scientific outputs. STC remains a dense method with an explicit
RAM budget. Version: 7.1.2.

## Prediction target (7.1.7)

**Restore automatic residual trend** is enabled by default, preserving previous
results. The observation reduction is always `y - custom - automaticTrend`.
The covariance/spline fit uses exactly that reduced signal in both modes.

- Enabled: `prediction = custom + automaticTrend + residual`.
- Disabled: `prediction = custom + residual`. This interprets the automatic
  field as a systematic observation contribution excluded from the desired
  deformation, a scientific assumption selected by the user.

The covariance shape alone cannot identify a physical systematic contribution.
The original custom model is always restored. Fixed custom parameters retain
their values. With exclusion, zero only the automatic field's target basis:
keep its observation design, fitting derivative and correlation with the
residual estimate in uncertainty propagation. Assigned parameter errors follow
the same target. Output uncertainty is conditional on this interpretation,
selected field/covariance, regularisation and outlier mask; it excludes model
selection, regularisation bias and physical model discrepancy.

Component CSVs retain `TemporalTrend_mm` for the full removed contribution and
add `RestoredTrend_mm` for the contribution actually applied. Thus
`Total_mm = CustomModel_mm + RestoredTrend_mm + Residual_mm`. JSON/MAT record
the restore flag and target equation. References remain unchanged.

The temporary temporal-example button has been removed. The normal expression,
parameter table and preview remain. Residual-mean/full assigned-covariance
editors remain hidden; existing profile settings are preserved and indicated.
No residual diagnostic figures are generated. Numeric checks remain available.

## Consolidated interface (7.1.8)

At the top of **Custom deterministic model**, choose **Enter expression** or
**Load MATLAB function (.m)**. Only the active source is displayed and evaluated.
Switching preserves the expression and file path, but an inactive file never
overrides an expression. Profiles created before this selector infer their
source from the previously saved function path. Parameters, coordinates and the
restore-trend choice are shared by both modes.

The full-width preview evaluates current fixed/initial values at illustrative
positions; it does not estimate parameters or run the modelling workflow.
Coordinate convention is followed by X/Y origins on the same row.

The AOI page shows the shapefile boundary below its path using the same reader
and CRS transformation as processing. File mode permits navigation and fitting
the extent, with drawing controls hidden. Drawing mode retains polygon editing.

The 1D/2D covariance sections show model, bin with optional manual value, and
suggested/maximum correlation length in that order. The secondary shared
refinement section bounds refinement around supplied suggestions; it has no
effect when both suggested lengths are empty. Output controls are grouped into
PS time-series export and query-point interpolation. These interface changes do
not modify the validated numerical modelling or uncertainty calculations.

## Excel report dependencies and repair (7.1.10)

PNG embedding requires both openpyxl and Pillow. Installation and the modelling
preflight verify actual PNG insertion into an XLSX using the report interpreter;
missing packages are prepared before the numerical model starts. An image or
formatting failure is reported as an incomplete report, with saved numerical
results retained.

For an existing report, close it in Excel and select **Outputs → Repair saved
Excel report**, then choose the run workbook. Keep its `figures` folder beside
it. This repairs the report without rerunning preprocessing, PSI or modelling.
A backup is kept under `files/report_backups`; the existing workbook is replaced
only after all repair steps succeed. Repeating repair replaces the PHASE figures
at their own anchors, preserving other images and numerical cell values.
