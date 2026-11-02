# PHASE 7.1 — Custom deterministic models

Automatic is the default. The four original automatic scientific functions
are unchanged. Each separate `STmodel_*_custom.m` entry copies its original
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

## Additional outputs

- `custom_model.json`: model, parameters/covariance, units, origins, references,
  native residual settings and uncertainty scope.
- `custom_model_results.mat`: retained/original observations, components,
  coordinates and variance decomposition.
- `custom_model_components.csv`: component values and signed covariance terms.
- `custom_residual_mean*.csv/json`: residual-mean diagnostics.
- `custom_model_source/`: copied custom MATLAB function, when selected.
- `figures/custom_covariance_time.*` and `custom_covariance_space.*`:
  PNG and editable FIG plots of empirical values, contributing pair counts,
  marginal fits and the common-amplitude kernel actually used for prediction.
  Matching CSVs retain the plotted empirical/kernel values and pair counts.
- `figures/custom_residual_mean.*`: mean/RMS by acquisition and mean by PS.
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
and mean-by-PS diagnostics separately. No additional detrending is performed.

The normal workbook and spatial export pipeline also receive custom results.

## Verification

```matlab
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
