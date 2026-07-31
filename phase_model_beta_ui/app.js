"use strict";

const ModelUI = {
  bridge: null,
  state: null,
  schema: { groups: [], fields: [] },
  config: {},
  active: "project",
  logs: [],
  localDirty: false,
  liveRaw: "",
  liveLines: [],
  liveTimer: null,
  clockTimer: null,
  clockReceivedAt: Date.now(),
  aoiMap: null,
};

function setup(htmlComponent) {
  ModelUI.bridge = htmlComponent;
  htmlComponent.addEventListener("DataChanged", () => receiveState(htmlComponent.Data));
  htmlComponent.addEventListener("PhaseLog", event => appendLog(event.Data));
  htmlComponent.addEventListener("PhaseStatus", event => receiveStatus(event.Data));
  htmlComponent.addEventListener("MapTilesReady", event => handleMapTilesReady(event.Data));
  document.addEventListener("DOMContentLoaded", wireStaticControls, { once: true });
  if (document.readyState !== "loading") wireStaticControls();
  htmlComponent.sendEventToMATLAB("Ready", {});
  if (htmlComponent.Data) receiveState(htmlComponent.Data);
}

function wireStaticControls() {
  if (wireStaticControls.done) return;
  wireStaticControls.done = true;
  listen("root-button", "click", () => send("OpenRoot", {}));
  listen("load-button", "click", () => send("Load", {}));
  listen("save-button", "click", save);
  listen("banner-save", "click", save);
  listen("start-button", "click", start);
  listen("stop-button", "click", () => send("Stop", {}));
  listen("clear-log", "click", () => {
    ModelUI.logs = [];
    ModelUI.liveLines = [];
    renderLogs();
    send("ClearLog", {});
  });
  listen("previous-section", "click", () => moveSection(-1));
  listen("next-section", "click", () => moveSection(1));
  ModelUI.aoiMap = new PhaseMap(byId("aoi-map"), {
    onPolygonChanged: handleMapPolygon,
    onTilesRequested: requests => {
      if (requests?.length) send("MapTilesRequested", { requests });
    },
    onDrawingChanged: drawing => {
      byId("map-draw").classList.toggle("hidden", drawing);
      byId("map-finish").classList.toggle("hidden", !drawing);
    },
  });
  listen("map-draw", "click", () => ModelUI.aoiMap.startDrawing());
  listen("map-finish", "click", () => ModelUI.aoiMap.finishDrawing());
  listen("map-fit", "click", () => {
    if (ModelUI.aoiMap.polygon.length) ModelUI.aoiMap.fitBounds(ModelUI.aoiMap.polygon);
  });
  listen("map-rectangle", "click", useCoordinateRectangle);
  if (!ModelUI.clockTimer) {
    ModelUI.clockTimer = window.setInterval(() => {
      if (ModelUI.state?.running) renderProgress();
    }, 1000);
  }
}

function receiveState(state) {
  if (!state || state.kind !== "state") return;
  const wasRunning = Boolean(ModelUI.state?.running);
  ModelUI.state = state;
  ModelUI.schema = state.schema || { groups: [], fields: [] };
  ModelUI.config = { ...(state.config || {}) };
  ModelUI.localDirty = Boolean(state.dirty);
  ModelUI.logs = asArray(state.logs);
  ModelUI.clockReceivedAt = Date.now();
  if (!wasRunning && state.running) {
    ModelUI.liveRaw = "";
    ModelUI.liveLines = [];
  }
  if (!visibleGroups().some(group => group.id === ModelUI.active) && ModelUI.active !== "run") {
    ModelUI.active = "project";
  }
  renderAll();
  syncLiveLogPolling();
}

function receiveStatus(status) {
  if (!status || !ModelUI.state) return;
  ModelUI.state = { ...ModelUI.state, ...status };
  ModelUI.localDirty = Boolean(status.dirty);
  ModelUI.clockReceivedAt = Date.now();
  renderStatus();
  renderProgress();
  syncLiveLogPolling();
}

function renderAll() {
  renderNavigation();
  renderPage();
  renderStatus();
  renderProgress();
  renderLogs();
  renderMap();
  byId("root-button").textContent = ModelUI.state?.rootDir || "PHASE project";
  byId("version").textContent = `PHASE Model ${ModelUI.state?.version || "Beta"}`;
}

function visibleGroups() {
  const dimension = ModelUI.config.projDim || "1D";
  const temporal = ["temporal", "temporal&NNI"].includes(ModelUI.config.procType);
  return asArray(ModelUI.schema.groups).filter(group => {
    if (group.id === "observations" || group.id === "splines") return temporal;
    if (group.id === "spatial1d") return !temporal && dimension === "1D";
    if (group.id === "spatial2d") return !temporal && dimension === "2D";
    return true;
  });
}

function renderNavigation() {
  const nav = byId("navigation");
  nav.innerHTML = `<div class="nav-label">Configuration</div>`;
  visibleGroups().forEach((group, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `nav-item${ModelUI.active === group.id ? " active" : ""}`;
    button.innerHTML = `<span class="nav-index">${index + 1}</span><span class="nav-title">${escapeHtml(group.title)}</span>`;
    button.addEventListener("click", () => switchPage(group.id));
    nav.appendChild(button);
  });
  nav.insertAdjacentHTML("beforeend", `<div class="nav-label">Tools</div>`);
  const monitor = document.createElement("button");
  monitor.type = "button";
  monitor.className = `nav-item${ModelUI.active === "run" ? " active" : ""}`;
  monitor.innerHTML = `<span class="nav-index">⌁</span><span class="nav-title">Run monitor</span>`;
  monitor.addEventListener("click", () => switchPage("run"));
  nav.appendChild(monitor);
}

function switchPage(id) {
  collectVisibleForm(false);
  ModelUI.active = id;
  renderNavigation();
  renderPage();
}

function renderPage() {
  byId("form-page").classList.toggle("active", ModelUI.active !== "run");
  byId("run-page").classList.toggle("active", ModelUI.active === "run");
  if (ModelUI.active === "run") {
    renderProgress();
    renderLogs();
    return;
  }
  const group = visibleGroups().find(item => item.id === ModelUI.active) || visibleGroups()[0];
  if (!group) return;
  ModelUI.active = group.id;
  byId("aoi-map-panel").classList.toggle("hidden", group.id !== "aoi");
  byId("section-kicker").textContent = group.kicker || "Configuration";
  byId("section-title").textContent = group.title || "";
  byId("section-subtitle").textContent = group.subtitle || "";
  const fields = asArray(ModelUI.schema.fields).filter(item => item.group === group.id);
  const relevant = fields.filter(isRelevant);
  byId("section-meta").textContent = `${relevant.length} control${relevant.length === 1 ? "" : "s"}`;
  const grid = byId("form-grid");
  grid.innerHTML = "";
  fields.forEach(meta => {
    const card = createField(meta);
    if (!isRelevant(meta)) card.classList.add("inactive");
    grid.appendChild(card);
    if (meta.id === "flag_AOIbb" && ModelUI.config.flag_AOIbb) {
      grid.appendChild(createFullPsExtentAction());
    }
  });
  renderContextNote(group.id);
  const groups = visibleGroups();
  const index = groups.findIndex(item => item.id === group.id);
  byId("previous-section").disabled = index <= 0;
  byId("next-section").disabled = index < 0 || index >= groups.length - 1;
  if (group.id === "aoi") window.setTimeout(() => ModelUI.aoiMap?.queueRender(),30);
}

function renderContextNote(groupId) {
  const note = byId("context-note");
  let message = "";
  if (groupId === "processing") {
    message = `${processingLabel()} · ${ModelUI.config.projDim || "1D"} project`;
  } else if (groupId === "spatial1d" || groupId === "spatial2d") {
    message = ModelUI.config.procType === "spatialDET"
      ? "Showing deterministic spline controls for the selected project dimension."
      : "Showing stochastic covariance controls for the selected project dimension.";
  }
  note.textContent = message;
  note.classList.toggle("hidden", !message);
}

function createField(meta) {
  const card = document.createElement("div");
  card.className = `field-card${meta.type === "path" ? " wide" : ""}`;
  if (meta.id === "flag_AOIbb") card.classList.add("wide");
  const role = pairRole(meta.id);
  if (role) card.classList.add(`pair-${role}`);
  const covarianceRole = covarianceModelRole(meta.id);
  if (covarianceRole) card.classList.add("covariance-model", covarianceRole);
  if (["lonMinAOI", "latMinAOI"].includes(meta.id)) card.classList.add("bounds-left");
  if (["lonMaxAOI", "latMaxAOI"].includes(meta.id)) card.classList.add("bounds-right");
  card.dataset.field = meta.id;
  const header = document.createElement("div");
  header.className = "field-label-row";
  header.innerHTML = `<label class="field-label" for="field-${escapeAttribute(meta.id)}">${escapeHtml(meta.label)}</label>` +
    `<span class="field-unit">${escapeHtml(meta.unit || "")}</span>`;
  card.appendChild(header);
  const row = document.createElement("div");
  row.className = "control-row";
  row.appendChild(createControl(meta));
  if (meta.browse) {
    const browse = document.createElement("button");
    browse.type = "button";
    browse.className = "browse-button";
    browse.textContent = "Browse…";
    browse.disabled = Boolean(ModelUI.state?.running);
    browse.addEventListener("click", () => {
      collectVisibleForm(false);
      send("Browse", { field: meta.id, config: { ...ModelUI.config } });
    });
    row.appendChild(browse);
  }
  card.appendChild(row);
  const description = document.createElement("div");
  description.className = "field-description";
  description.textContent = meta.description || "";
  card.appendChild(description);
  return card;
}

function createFullPsExtentAction() {
  const card = document.createElement("div");
  card.className = "field-card wide aoi-extent-card";
  const copy = document.createElement("div");
  copy.innerHTML = `<strong>Use the complete PS area</strong><span>Estimate a padded bounding box from every valid PS coordinate in the selected displacement file.</span>`;
  const button = document.createElement("button");
  button.type = "button";
  button.className = "button secondary";
  button.textContent = "Select full PS extent";
  button.disabled = Boolean(ModelUI.state?.running) || !textValue(ModelUI.config.filepathIN);
  button.addEventListener("click", () => {
    collectVisibleForm(false);
    send("EstimatePsBounds", payload());
  });
  card.append(copy, button);
  return card;
}

function createControl(meta) {
  const value = ModelUI.config[meta.id];
  if (meta.type === "boolean") {
    const wrap = document.createElement("div");
    wrap.className = "toggle-wrap";
    const toggle = document.createElement("button");
    toggle.type = "button";
    toggle.id = `field-${meta.id}`;
    toggle.className = `toggle${value ? " on" : ""}`;
    toggle.dataset.configField = meta.id;
    toggle.dataset.controlType = "boolean";
    toggle.setAttribute("role", "switch");
    toggle.setAttribute("aria-checked", String(Boolean(value)));
    toggle.disabled = Boolean(ModelUI.state?.running);
    const state = document.createElement("span");
    state.className = "toggle-state";
    state.textContent = value ? "Enabled" : "Disabled";
    toggle.addEventListener("click", () => {
      const next = toggle.getAttribute("aria-checked") !== "true";
      toggle.setAttribute("aria-checked", String(next));
      toggle.classList.toggle("on", next);
      state.textContent = next ? "Enabled" : "Disabled";
      ModelUI.config[meta.id] = next;
      localChange();
      if (dependencyDriver(meta.id)) renderAll();
    });
    wrap.append(toggle, state);
    return wrap;
  }

  let control;
  if (meta.type === "select") {
    control = document.createElement("select");
    const options = asArray(meta.options);
    const labels = asArray(meta.optionLabels);
    options.forEach((option, index) => {
      const node = document.createElement("option");
      node.value = textValue(option);
      node.textContent = textValue(labels[index] ?? option);
      control.appendChild(node);
    });
    control.value = textValue(value);
  } else {
    control = document.createElement("input");
    if (meta.type === "number") {
      control.type = "text";
      control.inputMode = "decimal";
      control.value = numberValue(value);
    } else if (meta.type === "date") {
      control.type = "date";
      control.value = dateInputValue(value);
    } else {
      control.type = "text";
      control.value = textValue(value);
      control.spellcheck = false;
    }
  }
  control.id = `field-${meta.id}`;
  control.className = "control";
  control.dataset.configField = meta.id;
  control.dataset.controlType = meta.type;
  control.disabled = Boolean(ModelUI.state?.running);
  control.addEventListener("change", () => {
    collectControl(control);
    localChange();
    if (dependencyDriver(meta.id)) renderAll();
  });
  control.addEventListener("input", () => {
    collectControl(control);
    localChange(false);
  });
  return control;
}

function dependencyDriver(id) {
  return [
    "flag_AOIbb", "projDim", "procType", "varNoise_method", "coll_proc",
    "num_spl_method", "lambda_method", "detrendMethod", "flag_tsExtr",
    "min_period_days_method", "min_coll_snr_method",
    "min_coll_corr_samples_method", "spline_min_knot_intervals_method",
    "spline_max_fraction_method", "dtCov_method_STC1D",
    "dsCov_method_STC1D", "dtCov_method_STC2D", "dsCov_method_STC2D",
    "varNoise_DET1D", "num_spl_method_DET1D", "lambda_method_DET1D",
    "varNoise_DET2D", "num_spl_method_DET2D", "lambda_method_DET2D",
  ].includes(id);
}

function isRelevant(meta) {
  const c = ModelUI.config;
  const spatial = ["spatialDET", "spatialSTC"].includes(c.procType);
  const needsInterpolationGrid = c.procType !== "temporal";
  const deterministic = c.procType === "spatialDET";
  const stochastic = c.procType === "spatialSTC";
  const rules = {
    filepathAOI: !c.flag_AOIbb,
    lonMinAOI: Boolean(c.flag_AOIbb), lonMaxAOI: Boolean(c.flag_AOIbb),
    latMinAOI: Boolean(c.flag_AOIbb), latMaxAOI: Boolean(c.flag_AOIbb),
    minMonths: spatial,
    cline_resolution: needsInterpolationGrid && c.projDim === "1D",
    grid_resolution: needsInterpolationGrid && c.projDim === "2D",
    step_t_ST: spatial,
    varNoise_manual: c.varNoise_method === "manual",
    coherence_dir: c.varNoise_method === "coherence",
    constellation: c.varNoise_method === "coherence",
    num_looks: c.varNoise_method === "coherence",
    coll_step_est: c.coll_proc === "prediction",
    min_coll_snr: c.min_coll_snr_method === "manual",
    min_coll_corr_samples: c.min_coll_corr_samples_method === "manual",
    spline_method: c.num_spl_method === "auto",
    num_spl_manual: c.num_spl_method === "manual",
    lambda_manual: c.lambda_method === "manual",
    min_period_days: c.min_period_days_method === "manual",
    spline_min_knot_intervals: c.spline_min_knot_intervals_method === "manual",
    spline_max_fraction: c.spline_max_fraction_method === "manual",
    detrendMethod: spatial, polyDegreeST: spatial,
    useInclinedMeansST: spatial && c.detrendMethod === "residualAtm",
    filepath_EXTR: Boolean(c.flag_tsExtr),
    dtCov_STC1D: stochastic && c.dtCov_method_STC1D === "manual",
    dsCov_STC1D: stochastic && c.dsCov_method_STC1D === "manual",
    dtCov_STC2D: stochastic && c.dtCov_method_STC2D === "manual",
    dsCov_STC2D: stochastic && c.dsCov_method_STC2D === "manual",
    varNoise_manual_DET1D: deterministic && c.varNoise_DET1D === "manual",
    spline_method_DET1D: deterministic && c.num_spl_method_DET1D === "auto",
    num_spl_row_manual_DET1D: deterministic && c.num_spl_method_DET1D === "manual",
    num_spl_col_manual_DET1D: deterministic && c.num_spl_method_DET1D === "manual",
    lambda_manual_DET1D: deterministic && c.lambda_method_DET1D === "manual",
    varNoise_manual_DET2D: deterministic && c.varNoise_DET2D === "manual",
    spline_method_DET2D: deterministic && c.num_spl_method_DET2D === "auto",
    num_spl_row_manual_DET2D: deterministic && c.num_spl_method_DET2D === "manual",
    num_spl_col_manual_DET2D: deterministic && c.num_spl_method_DET2D === "manual",
    num_spl_t_manual_DET2D: deterministic && c.num_spl_method_DET2D === "manual",
    lambda_manual_DET2D: deterministic && c.lambda_method_DET2D === "manual",
  };
  if (Object.prototype.hasOwnProperty.call(rules, meta.id)) return rules[meta.id];
  if (meta.group === "spatial1d") {
    if (meta.id.includes("_STC1D") || meta.id.includes("CovModel_STC1D")) return stochastic;
    return deterministic;
  }
  if (meta.group === "spatial2d") {
    if (meta.id.includes("_STC2D") || meta.id.includes("CovModel_STC2D")) return stochastic;
    return deterministic;
  }
  return true;
}

function pairRole(id) {
  const methods = new Set([
    "varNoise_method", "coll_proc", "num_spl_method", "lambda_method",
    "min_period_days_method", "min_coll_snr_method",
    "min_coll_corr_samples_method", "spline_min_knot_intervals_method",
    "spline_max_fraction_method", "dtCov_method_STC1D",
    "dsCov_method_STC1D", "dtCov_method_STC2D", "dsCov_method_STC2D",
    "varNoise_DET1D", "num_spl_method_DET1D", "lambda_method_DET1D",
    "varNoise_DET2D", "num_spl_method_DET2D", "lambda_method_DET2D",
  ]);
  const values = new Set([
    "varNoise_manual", "coherence_dir", "coll_step_est",
    "spline_method", "num_spl_manual", "lambda_manual",
    "min_period_days", "min_coll_snr", "min_coll_corr_samples",
    "spline_min_knot_intervals", "spline_max_fraction",
    "dtCov_STC1D", "dsCov_STC1D", "dtCov_STC2D", "dsCov_STC2D",
    "varNoise_manual_DET1D", "spline_method_DET1D",
    "num_spl_row_manual_DET1D", "num_spl_col_manual_DET1D",
    "lambda_manual_DET1D",
    "varNoise_manual_DET2D", "spline_method_DET2D",
    "num_spl_row_manual_DET2D", "num_spl_col_manual_DET2D",
    "num_spl_t_manual_DET2D", "lambda_manual_DET2D",
  ]);
  if (methods.has(id)) return "method";
  if (values.has(id)) return "value";
  return "";
}

function covarianceModelRole(id) {
  if (["tCovModel_STC1D", "tCovModel_STC2D"].includes(id)) return "covariance-temporal";
  if (["sCovModel_STC1D", "sCovModel_STC2D"].includes(id)) return "covariance-spatial";
  return "";
}

function collectVisibleForm(notify = true) {
  document.querySelectorAll("[data-config-field]").forEach(collectControl);
  if (notify) localChange();
}

function collectControl(control) {
  const id = control.dataset.configField;
  const type = control.dataset.controlType;
  if (!id) return;
  if (type === "boolean") {
    ModelUI.config[id] = control.getAttribute("aria-checked") === "true";
  } else if (type === "number") {
    ModelUI.config[id] = control.value.trim() === "" ? null : Number(control.value);
  } else {
    ModelUI.config[id] = control.value;
  }
  if (["lonMinAOI","lonMaxAOI","latMinAOI","latMaxAOI"].includes(id)) {
    ModelUI.config.aoi_polygon_lonlat = [];
  }
}

function localChange(notify = true) {
  ModelUI.localDirty = true;
  if (ModelUI.state) {
    ModelUI.state.dirty = true;
    if (!ModelUI.state.running) {
      ModelUI.state.status = "idle";
      ModelUI.state.statusDetail = "Unsaved changes";
    }
  }
  renderStatus();
  if (notify) send("Changed", payload());
}

function save() {
  collectVisibleForm(false);
  send("Save", payload());
}

function start() {
  collectVisibleForm(false);
  ModelUI.active = "run";
  renderNavigation();
  renderPage();
  send("Start", payload());
}

function payload() { return { config: { ...ModelUI.config } }; }

function renderStatus() {
  const state = ModelUI.state || {};
  const status = state.status || "idle";
  byId("status-dot").className = `status-dot ${status}`;
  byId("status-title").textContent = status === "idle" && (state.dirty || ModelUI.localDirty)
    ? "Unsaved" : status;
  byId("status-detail").textContent = state.statusDetail || "Waiting for MATLAB";
  byId("dirty-banner").classList.toggle("hidden", !(state.dirty || ModelUI.localDirty));
  const running = Boolean(state.running);
  byId("load-button").disabled = running;
  byId("save-button").disabled = running;
  byId("start-button").disabled = running || state.dirty || ModelUI.localDirty || status === "error";
  byId("stop-button").disabled = !running;
}

function renderMap() {
  if (!ModelUI.aoiMap || !ModelUI.state?.map) return;
  const polygon = matrix(ModelUI.state.map.polygon?.length
    ? ModelUI.state.map.polygon
    : ModelUI.config.aoi_polygon_lonlat);
  ModelUI.aoiMap.setData({
    coastlines: asArray(ModelUI.state.map.coastlines),
    footprints: [],
    polygon,
  });
  updateMapMeta(polygon);
}

function handleMapTilesReady(data) {
  ModelUI.aoiMap?.retryTiles(asArray(data?.keys));
}

function handleMapPolygon(polygon) {
  const closed = closePolygon(polygon);
  if (closed.length < 4) return;
  const lons = closed.map(point => Number(point[0]));
  const lats = closed.map(point => Number(point[1]));
  const bbox = {
    minLon: Math.min(...lons), maxLon: Math.max(...lons),
    minLat: Math.min(...lats), maxLat: Math.max(...lats),
  };
  ModelUI.config.flag_AOIbb = true;
  ModelUI.config.aoi_polygon_lonlat = closed;
  ModelUI.config.lonMinAOI = bbox.minLon;
  ModelUI.config.lonMaxAOI = bbox.maxLon;
  ModelUI.config.latMinAOI = bbox.minLat;
  ModelUI.config.latMaxAOI = bbox.maxLat;
  localChange(false);
  updateMapMeta(closed);
  send("MapAoiChanged", { polygon: closed, bbox });
}

function useCoordinateRectangle() {
  collectVisibleForm(false);
  const values = ["lonMinAOI","lonMaxAOI","latMinAOI","latMaxAOI"]
    .map(name => Number(ModelUI.config[name]));
  if (!values.every(Number.isFinite) || values[0] >= values[1] || values[2] >= values[3]) return;
  const polygon = closePolygon([
    [values[0],values[2]],[values[1],values[2]],
    [values[1],values[3]],[values[0],values[3]],
  ]);
  ModelUI.aoiMap.setPolygon(polygon,true);
  handleMapPolygon(polygon);
}

function updateMapMeta(polygon) {
  const points = matrix(polygon);
  if (!points.length) {
    byId("map-bbox").textContent = "No AOI drawn";
    byId("map-vertices").textContent = "—";
    return;
  }
  const lons = points.map(point => point[0]), lats = points.map(point => point[1]);
  byId("map-bbox").textContent = `${Math.min(...lons).toFixed(6)}, ${Math.min(...lats).toFixed(6)} → ${Math.max(...lons).toFixed(6)}, ${Math.max(...lats).toFixed(6)}`;
  byId("map-vertices").textContent = `${Math.max(0,points.length-1)} vertices`;
}

function renderProgress() {
  const progress = ModelUI.state?.progress || { percentage: 0, phase: "Ready", elapsedSeconds: 0 };
  const percentage = Math.max(0, Math.min(100, Number(progress.percentage) || 0));
  let elapsed = Number(progress.elapsedSeconds) || 0;
  if (ModelUI.state?.running) elapsed += (Date.now() - ModelUI.clockReceivedAt) / 1000;
  byId("sidebar-progress").style.width = `${percentage}%`;
  byId("sidebar-progress-label").textContent = progress.phase || "Model ready";
  byId("run-progress").style.width = `${percentage}%`;
  byId("run-percentage").textContent = `${Math.round(percentage)}%`;
  byId("run-phase").textContent = progress.phase || "Ready";
  byId("run-elapsed").textContent = formatDuration(elapsed);
  byId("run-model-summary").textContent = `${processingLabel()} · ${ModelUI.config.projDim || "1D"}`;
}

function processingLabel() {
  const labels = {
    temporal: "Temporal",
    "temporal&NNI": "Temporal + NNI",
    spatialDET: "Spatio-temporal deterministic",
    spatialSTC: "Spatio-temporal stochastic",
  };
  return labels[ModelUI.config.procType] || textValue(ModelUI.config.procType);
}

function appendLog(entry) {
  if (!entry) return;
  ModelUI.logs.push(entry);
  if (ModelUI.logs.length > 800) ModelUI.logs = ModelUI.logs.slice(-800);
  if (!ModelUI.liveLines.length) renderLogs();
}

function renderLogs() {
  const consoleNode = byId("console");
  consoleNode.innerHTML = "";
  if (ModelUI.liveLines.length) {
    ModelUI.liveLines.forEach(message => consoleNode.appendChild(logElement({ time: "", message })));
  } else {
    ModelUI.logs.forEach(entry => consoleNode.appendChild(logElement(entry)));
  }
  const count = ModelUI.liveLines.length || ModelUI.logs.length;
  byId("log-count").textContent = `${count} message${count === 1 ? "" : "s"}`;
  consoleNode.scrollTop = consoleNode.scrollHeight;
}

function logElement(entry) {
  const line = document.createElement("div");
  line.className = "log-line";
  const time = typeof entry === "string" ? "" : textValue(entry.time);
  const message = typeof entry === "string" ? entry : textValue(entry.message);
  line.innerHTML = `<span class="log-time">${escapeHtml(time)}</span><span>${escapeHtml(message)}</span>`;
  return line;
}

function syncLiveLogPolling() {
  const url = ModelUI.state?.liveLogUrl;
  if (!url) return;
  if (ModelUI.state?.running && !ModelUI.liveTimer) {
    pollLiveLog();
    ModelUI.liveTimer = window.setInterval(pollLiveLog, 500);
  } else if (!ModelUI.state?.running && ModelUI.liveTimer) {
    window.clearInterval(ModelUI.liveTimer);
    ModelUI.liveTimer = null;
    pollLiveLog();
    window.setTimeout(pollLiveLog, 350);
  }
}

async function pollLiveLog() {
  const url = ModelUI.state?.liveLogUrl;
  if (!url) return;
  try {
    const response = await fetch(`${url}?v=${Date.now()}`, { cache: "no-store" });
    if (!response.ok) return;
    const raw = await response.text();
    if (raw === ModelUI.liveRaw) return;
    ModelUI.liveRaw = raw;
    ModelUI.liveLines = raw
      .replace(/\u001b\[[0-9;]*m/g, "")
      .split(/\r?\n/)
      .map(line => line.trimEnd())
      .filter(line => line.trim().length > 0);
    renderLogs();
  } catch (_) {
    // The diary file appears only after Start.
  }
}

function moveSection(delta) {
  collectVisibleForm(false);
  const groups = visibleGroups();
  const index = groups.findIndex(group => group.id === ModelUI.active);
  if (groups[index + delta]) switchPage(groups[index + delta].id);
}

function send(name, data) {
  if (ModelUI.bridge) ModelUI.bridge.sendEventToMATLAB(name, data || {});
}

function listen(id, eventName, callback) {
  const node = byId(id);
  if (node) node.addEventListener(eventName, callback);
}
function byId(id) { return document.getElementById(id); }
function asArray(value) { if (!value) return []; return Array.isArray(value) ? value : [value]; }
function matrix(value) {
  if (!Array.isArray(value)) return [];
  let rows = value;
  while (rows.length === 1 && Array.isArray(rows[0]) && Array.isArray(rows[0][0])) rows = rows[0];
  if (rows.length === 2 && rows[0]?.length > 2 && rows[0].length === rows[1]?.length) {
    rows = rows[0].map((lon,index) => [lon,rows[1][index]]);
  }
  return rows.filter(row => Array.isArray(row) && row.length >= 2)
    .map(row => [Number(row[0]),Number(row[1])])
    .filter(row => Number.isFinite(row[0]) && Number.isFinite(row[1]));
}
function closePolygon(value) {
  const points = matrix(value).map(point => [...point]);
  if (!points.length) return [];
  const first = points[0], last = points[points.length-1];
  if (first[0] !== last[0] || first[1] !== last[1]) points.push([...first]);
  return points;
}
function textValue(value) { return value === null || value === undefined ? "" : String(value); }
function numberValue(value) { return value === null || value === undefined || Number.isNaN(Number(value)) ? "" : String(value); }
function dateInputValue(value) {
  const raw = textValue(value);
  if (/^\d{4}-\d{2}-\d{2}$/.test(raw)) return raw;
  const compact = raw.replace(/[^0-9]/g, "");
  return compact.length === 8 ? `${compact.slice(0,4)}-${compact.slice(4,6)}-${compact.slice(6,8)}` : "";
}
function formatDuration(value) {
  const seconds = Math.max(0, Math.floor(Number(value) || 0));
  const hours = Math.floor(seconds / 3600);
  const minutes = Math.floor((seconds % 3600) / 60);
  const rest = seconds % 60;
  return hours > 0
    ? `${String(hours).padStart(2,"0")}:${String(minutes).padStart(2,"0")}:${String(rest).padStart(2,"0")}`
    : `${String(minutes).padStart(2,"0")}:${String(rest).padStart(2,"0")}`;
}
function escapeHtml(value) {
  const node = document.createElement("div");
  node.textContent = textValue(value);
  return node.innerHTML;
}
function escapeAttribute(value) { return escapeHtml(value).replace(/`/g, "&#96;"); }

window.addEventListener("error", event => {
  const banner = byId("error-banner");
  if (!banner) return;
  banner.textContent = `Interface error: ${event.message}`;
  banner.classList.remove("hidden");
});
