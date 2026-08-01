"use strict";

const PhaseUI = {
  bridge: null,
  state: null,
  config: {},
  schema: { groups: [], items: [] },
  active: "images",
  localDirty: false,
  logs: [],
  aoiMap: null,
  downloadMap: null,
  downloaderOpen: false,
  updateOpen: false,
  selectedAsfNames: new Set(),
  selectedUpdateNames: new Set(),
  runProgress: {
    percentage: 0, phase: "Ready", elapsedSeconds: 0,
    etaSeconds: NaN, startedAt: "", indeterminate: false,
  },
  runClock: null,
  asfResultSignature: "",
  updateResultSignature: "",
  asfSortKey: "date",
  asfSortDirection: "asc",
  updateSortKey: "date",
  updateSortDirection: "asc",
};

function setup(htmlComponent) {
  PhaseUI.bridge = htmlComponent;
  htmlComponent.addEventListener("DataChanged", () => receiveState(htmlComponent.Data));
  htmlComponent.addEventListener("PhaseLog", event => appendLog(event.Data));
  htmlComponent.addEventListener("PhaseRunProgress", event => receiveRunProgress(event.Data));
  htmlComponent.addEventListener("MapTilesReady", event => handleMapTilesReady(event.Data));
  document.addEventListener("DOMContentLoaded", wireStaticControls, { once: true });
  if (document.readyState !== "loading") wireStaticControls();
  htmlComponent.sendEventToMATLAB("Ready", {});
  if (htmlComponent.Data) receiveState(htmlComponent.Data);
}

function wireStaticControls() {
  if (wireStaticControls.done) return;
  wireStaticControls.done = true;
  listen("load-button", "click", () => send("Load", {}));
  listen("save-button", "click", save);
  listen("banner-save", "click", save);
  listen("start-button", "click", start);
  listen("stop-button", "click", () => send("Stop", {}));
  listen("workdir", "click", () => send("OpenWorkDir", {}));
  listen("clear-log", "click", () => { PhaseUI.logs = []; renderLogs(); });
  listen("previous-section", "click", () => moveSection(-1));
  listen("next-section", "click", () => moveSection(1));
  listen("open-downloader", "click", () => toggleWorkflow("downloader", true));
  listen("close-downloader", "click", () => toggleWorkflow("downloader", false));
  listen("open-update", "click", () => toggleWorkflow("update", true));
  listen("close-update", "click", () => toggleWorkflow("update", false));
  listen("import-images", "click", () => { collectVisibleForm(); send("ImportImages", payload()); });
  listen("open-slaves-folder", "click", () => { collectVisibleForm(); send("OpenSlavesFolder", payload()); });
  listen("refresh-slaves", "click", () => send("RefreshSlaves", {}));

  PhaseUI.aoiMap = new PhaseMap(byId("aoi-map"), {
    onPolygonChanged: handleMapPolygon,
    onTilesRequested: requestMapTiles,
    onDrawingChanged: drawing => {
      byId("map-draw").classList.toggle("hidden", drawing);
      byId("map-finish").classList.toggle("hidden", !drawing);
    },
  });
  listen("map-draw", "click", () => PhaseUI.aoiMap.startDrawing());
  listen("map-finish", "click", () => PhaseUI.aoiMap.finishDrawing());
  listen("map-clear", "click", resetMapAoi);
  listen("map-fit", "click", () => PhaseUI.aoiMap.fitToFootprints());
  listen("map-refresh", "click", () => send("RefreshMap", {}));

  PhaseUI.downloadMap = new PhaseMap(byId("download-map"), {
    onPolygonChanged: handleDownloadPolygon,
    onTilesRequested: requestMapTiles,
    onDrawingChanged: drawing => {
      byId("download-map-draw").classList.toggle("hidden", drawing);
      byId("download-map-finish").classList.toggle("hidden", !drawing);
    },
  });
  listen("download-map-draw", "click", () => PhaseUI.downloadMap.startDrawing());
  listen("download-map-finish", "click", () => PhaseUI.downloadMap.finishDrawing());
  listen("download-map-fit", "click", () => PhaseUI.downloadMap.fitToFootprints());
  listen("download-use-processing-aoi", "click", useProcessingAoiForDownload);
  listen("reset-asf-filters", "click", resetAsfFilters);
  listen("search-asf", "click", searchAsf);
  listen("earthdata-login", "click", loginEarthdata);
  listen("earthdata-password", "keydown", event => { if (event.key === "Enter") loginEarthdata(); });
  listen("earthdata-logout", "click", () => send("DownloadLogout", {}));
  listen("asf-select-all", "change", event => selectAllAsf(event.target.checked));
  listen("download-selected", "click", downloadSelectedAsf);
  listen("asf-sort-key", "change", event => {
    PhaseUI.asfSortKey = event.target.value;
    renderDownloader();
  });
  listen("asf-sort-direction", "click", () => {
    PhaseUI.asfSortDirection = PhaseUI.asfSortDirection === "asc" ? "desc" : "asc";
    renderDownloader();
  });
  listen("stop-download-transfer", "click", () => send("StopDownload", {}));
  listen("refresh-update", "click", () => send("RefreshUpdate", {}));
  listen("search-update", "click", searchUpdate);
  listen("update-select-all", "change", event => selectAllUpdate(event.target.checked));
  listen("download-update", "click", downloadSelectedUpdate);
  listen("update-sort-key", "change", event => {
    PhaseUI.updateSortKey = event.target.value;
    renderUpdate();
  });
  listen("update-sort-direction", "click", () => {
    PhaseUI.updateSortDirection = PhaseUI.updateSortDirection === "asc" ? "desc" : "asc";
    renderUpdate();
  });
  listen("stop-update-transfer", "click", () => send("StopDownload", {}));
  document.querySelectorAll(".filter-disclosure").forEach(disclosure => {
    disclosure.addEventListener("toggle", () => {
      window.setTimeout(() => PhaseUI.downloadMap?.queueRender(), 30);
    });
  });
  if (!PhaseUI.runClock) {
    PhaseUI.runClock = window.setInterval(tickRunClock, 1000);
  }
}

function requestMapTiles(requests) {
  if (requests?.length) send("MapTilesRequested", { requests });
}

function handleMapTilesReady(data) {
  const keys = asArray(data?.keys);
  PhaseUI.aoiMap?.retryTiles(keys);
  PhaseUI.downloadMap?.retryTiles(keys);
  if (Number(data?.failedCount || 0) > 0 && !keys.length) {
    toast("Satellite background unavailable; the PHASE vector map remains active.");
  }
}

function receiveState(state) {
  if (!state || state.kind !== "state") return;
  PhaseUI.state = state;
  PhaseUI.schema = state.schema || { groups: [], items: [] };
  PhaseUI.config = { ...(state.config || {}) };
  PhaseUI.localDirty = Boolean(state.dirty);
  PhaseUI.logs = asArray(state.logs);
  PhaseUI.runProgress = { ...PhaseUI.runProgress, ...(state.runProgress || {}) };
  if (!findGroup(PhaseUI.active) && PhaseUI.active !== "run") PhaseUI.active = "images";
  renderAll();
}

function renderAll() {
  renderNavigation();
  renderPage();
  renderStatus();
  renderLogs();
  renderSummary();
  renderRunProgress();
  renderMap();
  byId("workdir").textContent = PhaseUI.state?.workDir || "No preprocessing folder";
  byId("version").textContent = `PHASE Preprocessing ${PhaseUI.state?.version || "6.0.0"}`;
}

function renderNavigation() {
  const nav = byId("navigation");
  nav.innerHTML = `<div class="nav-label">Preprocessing</div>`;
  PhaseUI.schema.groups.forEach((group, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `nav-item${PhaseUI.active === group.id ? " active" : ""}`;
    button.innerHTML = `<span class="nav-index">${index + 1}</span><span class="nav-text">${escapeHtml(group.title)}</span>`;
    button.addEventListener("click", () => switchPage(group.id));
    nav.appendChild(button);
  });
  nav.insertAdjacentHTML("beforeend", `<div class="nav-label">Tools</div>`);
  nav.appendChild(toolNav("run", "⌁", "Run monitor"));
}

function toolNav(id, icon, title) {
  const button = document.createElement("button");
  button.type = "button";
  button.className = `nav-item${PhaseUI.active === id ? " active" : ""}`;
  button.innerHTML = `<span class="nav-index">${icon}</span><span class="nav-text">${title}</span>`;
  button.addEventListener("click", () => switchPage(id));
  return button;
}

function switchPage(id) {
  collectVisibleForm();
  PhaseUI.active = id;
  renderNavigation();
  renderPage();
}

function renderPage() {
  ["form-page", "run-page"].forEach(id => byId(id).classList.remove("active"));
  if (PhaseUI.active === "run") { byId("run-page").classList.add("active"); renderSummary(); return; }
  byId("form-page").classList.add("active");
  const group = findGroup(PhaseUI.active) || PhaseUI.schema.groups[0];
  if (!group) return;
  byId("section-title").textContent = group.title;
  byId("section-subtitle").textContent = group.subtitle || "";
  const code = constellationCode();
  const items = PhaseUI.schema.items.filter(item => item.group === group.id && (!item.appliesTo || item.appliesTo === "both" || item.appliesTo === code));
  byId("section-meta").textContent = `${items.length} controls`;
  const grid = byId("form-grid");
  grid.innerHTML = "";
  items.forEach(item => grid.appendChild(createField(item)));
  applyFieldDependencies();
  updateResumeWarning();
  byId("images-panel").classList.toggle("hidden", group.id !== "images");
  byId("aoi-map-panel").classList.toggle("hidden", group.id !== "aoi");
  if (group.id === "images") renderImagesPanel();
  if (group.id === "aoi") renderMap();
  const index = PhaseUI.schema.groups.findIndex(candidate => candidate.id === group.id);
  byId("previous-section").disabled = index <= 0;
  byId("next-section").disabled = index >= PhaseUI.schema.groups.length - 1;
}

function createField(item) {
  const wrapper = document.createElement("div");
  const pairedDemPath = ["dem_file", "dem_file_coreg"].includes(item.id);
  const fullWidth = ["first_step", "num_gcp"].includes(item.id)
    || (item.type === "path" && !pairedDemPath);
  wrapper.className = `field${fullWidth ? " full" : ""}`;
  wrapper.dataset.field = item.id;
  const heading = document.createElement("div");
  heading.className = "field-label-row";
  heading.innerHTML = `<label for="field-${item.id}">${escapeHtml(item.label)}</label>${item.unit ? `<span class="unit">${escapeHtml(item.unit)}</span>` : ""}`;
  wrapper.appendChild(heading);

  let control;
  if (item.type === "choice") {
    control = document.createElement("select");
    asArray(item.options).forEach(option => {
      const element = document.createElement("option");
      const label = String(option);
      element.value = item.id === "first_step" ? label.split(/[—-]/)[0].trim() : label;
      element.textContent = label;
      control.appendChild(element);
    });
    control.value = textValue(PhaseUI.config[item.id]);
  } else if (item.type === "toggle") {
    const switchWrap = document.createElement("div");
    switchWrap.className = "switch-wrap";
    switchWrap.innerHTML = `<label class="switch"><input id="field-${item.id}" type="checkbox"><span class="switch-slider"></span></label><span class="switch-state"></span>`;
    control = switchWrap.querySelector("input");
    control.checked = Boolean(PhaseUI.config[item.id]);
    updateSwitchState(switchWrap, control.checked);
    control.addEventListener("change", () => updateSwitchState(switchWrap, control.checked));
    wrapper.appendChild(switchWrap);
  } else {
    control = document.createElement("input");
    control.type = "text";
    control.value = textValue(PhaseUI.config[item.id]);
    control.spellcheck = false;
  }
  control.id = `field-${item.id}`;
  control.dataset.configField = item.id;
  control.addEventListener("input", markDirty);
  control.addEventListener("change", () => {
    markDirty();
    applyFieldDependencies();
    if (item.id === "constellation") renderPage();
  });

  if (item.type === "path") {
    const row = document.createElement("div");
    row.className = "path-control";
    row.appendChild(control);
    const browse = document.createElement("button");
    browse.type = "button"; browse.className = "browse-button"; browse.textContent = "Browse…";
    browse.addEventListener("click", () => { collectVisibleForm(); send("Browse", { ...payload(), field: item.id }); });
    row.appendChild(browse); wrapper.appendChild(row);
  } else if (item.type !== "toggle") {
    wrapper.appendChild(control);
  }
  const help = document.createElement("div");
  help.className = "field-help"; help.textContent = item.help || " "; wrapper.appendChild(help);
  return wrapper;
}

function updateSwitchState(wrapper, checked) {
  wrapper.querySelector(".switch-state").textContent = checked ? "Enabled" : "Disabled";
}

function collectVisibleForm() {
  document.querySelectorAll("[data-config-field]").forEach(control => {
    PhaseUI.config[control.dataset.configField] = control.type === "checkbox" ? control.checked : control.value;
  });
  applyAutomaticEpsg();
}

function markDirty() {
  collectVisibleForm();
  applyFieldDependencies();
  PhaseUI.localDirty = true;
  renderStatus(); renderNavigation(); renderSummary();
}

function applyFieldDependencies() {
  const resumeSlaves = Number(PhaseUI.config.first_step || 1) > 1;
  if (resumeSlaves) {
    PhaseUI.config.process_master = false;
    const processMaster = byId("field-process_master");
    if (processMaster) {
      processMaster.checked = false;
      const switchWrap = processMaster.closest(".switch-wrap");
      if (switchWrap) updateSwitchState(switchWrap,false);
    }
  }
  setFieldDisabled("master_date", Boolean(PhaseUI.config.auto_master));
  setFieldDisabled("process_master", resumeSlaves);
  setFieldDisabled("dem_file", PhaseUI.config.dem_name !== "External DEM");
  setFieldDisabled("dem_file_coreg", PhaseUI.config.dem_name_coreg !== "External DEM");
  setFieldDisabled("epsg_code", Boolean(PhaseUI.config.auto_epsg), true);
  updateResumeWarning();
}

function updateResumeWarning() {
  const warning = byId("resume-warning");
  if (!warning) return;
  const resumeSlaves = Number(PhaseUI.config.first_step || 1) > 1;
  const relevantPage = ["master","slaves"].includes(PhaseUI.active);
  warning.classList.toggle("hidden", !resumeSlaves || !relevantPage);
  if (resumeSlaves) {
    byId("resume-warning-step").textContent = textValue(PhaseUI.config.first_step);
  }
}

function setFieldDisabled(id, disabled, automatic = false) {
  const control = byId(`field-${id}`);
  if (!control) return;
  control.disabled = disabled;
  const wrapper = control.closest(".field");
  wrapper?.classList.toggle("inactive", disabled);
  wrapper?.classList.toggle("detected", disabled && automatic);
  const browse = wrapper?.querySelector(".browse-button");
  if (browse) browse.disabled = disabled;
}

function applyAutomaticEpsg() {
  if (!Boolean(PhaseUI.config.auto_epsg)) return;
  const code = estimateEpsg(
    Number(PhaseUI.config.lon_min), Number(PhaseUI.config.lat_min),
    Number(PhaseUI.config.lon_max), Number(PhaseUI.config.lat_max)
  );
  if (!code) return;
  PhaseUI.config.epsg_code = String(code);
  const control = byId("field-epsg_code");
  if (control) control.value = String(code);
}

function estimateEpsg(lonMin, latMin, lonMax, latMax) {
  if (![lonMin, latMin, lonMax, latMax].every(Number.isFinite)
      || lonMin >= lonMax || latMin >= latMax) return 0;
  const lon = (lonMin + lonMax) / 2;
  const lat = (latMin + latMax) / 2;
  if (lat >= 84) return 3413;
  if (lat <= -80) return 3031;
  let zone = Math.min(60, Math.max(1, Math.floor((lon + 180) / 6) + 1));
  if (lat >= 56 && lat < 64 && lon >= 3 && lon < 12) zone = 32;
  else if (lat >= 72 && lat < 84) {
    if (lon >= 0 && lon < 9) zone = 31;
    else if (lon < 21) zone = 33;
    else if (lon < 33) zone = 35;
    else if (lon < 42) zone = 37;
  }
  return (lat >= 0 ? 32600 : 32700) + zone;
}

function renderStatus() {
  const status = PhaseUI.state?.status || "idle";
  const detail = PhaseUI.state?.statusDetail || "Ready";
  const dirty = PhaseUI.localDirty;
  byId("status-dot").className = `status-dot ${status}`;
  byId("status-title").textContent = dirty && status !== "running" ? "Unsaved" : status;
  byId("status-detail").textContent = dirty && status !== "running" ? "Save before starting" : detail;
  byId("dirty-banner").classList.toggle("hidden", !dirty || PhaseUI.state?.running);
  byId("start-button").disabled = dirty || PhaseUI.state?.running;
  byId("save-button").disabled = PhaseUI.state?.running;
  byId("load-button").disabled = PhaseUI.state?.running;
  byId("stop-button").disabled = !PhaseUI.state?.running;
  const first = Number(PhaseUI.config.first_step || 1);
  const width = Math.max(0, Math.min(100, ((7 - first) / 6) * 100));
  byId("range-progress").style.marginLeft = `${100 - width}%`;
  byId("range-progress").style.width = `${width}%`;
  byId("range-label").textContent = `${PhaseUI.config.constellation || "Satellite"} · ${stepRangeLabel(first)}`;
}

function stepRangeLabel(first) {
  const descriptions = {
    1: "Prepare slaves", 2: "Split / subset", 3: "Coregister + IFG",
    4: "StaMPS export", 5: "Average intensity", 6: "Terrain correction",
  };
  return `${first} ${descriptions[first] || "Resume"} → 6 Terrain correction`;
}

function renderSummary() {
  byId("summary-constellation").textContent = PhaseUI.config.constellation || "—";
  byId("summary-first").textContent = textValue(PhaseUI.config.first_step) || "—";
  byId("summary-master").textContent = PhaseUI.config.auto_master ? "Automatic" : (PhaseUI.config.master_date || "Manual");
  byId("summary-epsg").textContent = PhaseUI.config.epsg_code ? `EPSG:${PhaseUI.config.epsg_code}` : "—";
}

function receiveRunProgress(progress) {
  if (!progress) return;
  PhaseUI.runProgress = { ...PhaseUI.runProgress, ...progress };
  renderRunProgress();
  if (PhaseUI.state?.running && byId("status-detail")) {
    byId("status-detail").textContent = PhaseUI.runProgress.phase || "Processing";
  }
}

function tickRunClock() {
  if (!PhaseUI.state?.running) return;
  PhaseUI.runProgress.elapsedSeconds = Number(PhaseUI.runProgress.elapsedSeconds || 0) + 1;
  const eta = Number(PhaseUI.runProgress.etaSeconds);
  if (Number.isFinite(eta) && eta >= 0) {
    PhaseUI.runProgress.etaSeconds = Math.max(0,eta - 1);
  }
  renderRunProgress();
}

function renderRunProgress() {
  const progress = PhaseUI.runProgress || {};
  const percentage = Math.max(0,Math.min(100,Number(progress.percentage) || 0));
  const bar = byId("processing-progress-bar");
  if (!bar) return;
  bar.style.width = `${percentage}%`;
  bar.classList.toggle("indeterminate", Boolean(progress.indeterminate));
  byId("processing-progress-percent").textContent = `${Math.round(percentage)}%`;
  byId("processing-progress-phase").textContent = progress.phase || "Ready";
  byId("processing-elapsed").textContent = formatDuration(progress.elapsedSeconds);
  const eta = Number(progress.etaSeconds);
  byId("processing-eta").textContent = Number.isFinite(eta) && eta >= 0
    ? formatDuration(eta)
    : (PhaseUI.state?.running ? "Estimating…" : "—");
}

function formatDuration(rawSeconds) {
  const total = Math.max(0,Math.round(Number(rawSeconds) || 0));
  const hours = Math.floor(total / 3600);
  const minutes = Math.floor((total % 3600) / 60);
  const seconds = total % 60;
  if (hours) return `${hours}h ${String(minutes).padStart(2,"0")}m ${String(seconds).padStart(2,"0")}s`;
  return `${minutes}m ${String(seconds).padStart(2,"0")}s`;
}

function renderImagesPanel() {
  const isSEN = constellationCode() === "SEN";
  byId("sentinel-tools").classList.toggle("hidden", !isSEN);
  if (!isSEN) {
    PhaseUI.downloaderOpen = false;
    PhaseUI.updateOpen = false;
  }
  byId("import-description").textContent = isSEN
    ? "Move or copy Sentinel-1 ZIP products into the project slaves folder."
    : "Copy COSMO-SkyMed HDF5 products into the project slaves folder.";
  byId("downloader-panel").classList.toggle("hidden", !isSEN || !PhaseUI.downloaderOpen);
  byId("update-panel").classList.toggle("hidden", !isSEN || !PhaseUI.updateOpen);
  renderSlaveInventory();
  if (isSEN) {
    renderDownloader();
    renderUpdate();
    renderTransfer();
  }
}

function renderMap() {
  if (!PhaseUI.aoiMap || !PhaseUI.state?.map) return;
  const mapData = PhaseUI.state.map;
  const footprintLabel = constellationCode() === "SEN" ? "Local Sentinel-1" : "Local COSMO-SkyMed";
  const footprints = asArray(mapData.footprints).filter(item => String(item?.source || "") === footprintLabel);
  PhaseUI.aoiMap.setData({
    coastlines: asArray(mapData.coastlines),
    footprints,
    polygon: matrix(mapData.polygon),
  });
  const polygon = matrix(mapData.polygon);
  updateMapMeta(polygon, footprints.length);
}

function handleMapPolygon(polygon) {
  if (!polygon || polygon.length < 3) return;
  const closed = closePolygon(polygon);
  const lons = closed.map(point => Number(point[0]));
  const lats = closed.map(point => Number(point[1]));
  const bbox = {
    minLon: Math.min(...lons), maxLon: Math.max(...lons),
    minLat: Math.min(...lats), maxLat: Math.max(...lats),
  };
  PhaseUI.config.lon_min = formatCoordinate(bbox.minLon);
  PhaseUI.config.lon_max = formatCoordinate(bbox.maxLon);
  PhaseUI.config.lat_min = formatCoordinate(bbox.minLat);
  PhaseUI.config.lat_max = formatCoordinate(bbox.maxLat);
  applyAutomaticEpsg();
  PhaseUI.localDirty = true;
  updateMapMeta(closed, PhaseUI.aoiMap.footprints.length);
  renderStatus(); renderSummary(); renderNavigation();
  send("MapAoiChanged", { polygon: closed, bbox });
}

function resetMapAoi() {
  const polygon = bboxPolygonFromConfig();
  PhaseUI.aoiMap.setPolygon(polygon, true);
  handleMapPolygon(polygon);
}

function updateMapMeta(polygon, footprintCount) {
  if (polygon && polygon.length) {
    const lons = polygon.map(point => Number(point[0]));
    const lats = polygon.map(point => Number(point[1]));
    byId("map-bbox").textContent = `${formatCoordinate(Math.min(...lons))}, ${formatCoordinate(Math.min(...lats))} → ${formatCoordinate(Math.max(...lons))}, ${formatCoordinate(Math.max(...lats))}`;
  } else {
    byId("map-bbox").textContent = "No AOI drawn";
  }
  byId("map-footprint-count").textContent = footprintCount
    ? `${footprintCount} footprint${footprintCount === 1 ? "" : "s"} from slaves`
    : "No footprint in slaves";
}

function toggleWorkflow(kind, open) {
  if (kind === "downloader") {
    PhaseUI.downloaderOpen = open;
    if (open) PhaseUI.updateOpen = false;
  } else {
    PhaseUI.updateOpen = open;
    if (open) PhaseUI.downloaderOpen = false;
  }
  renderImagesPanel();
  if (open) window.setTimeout(() => {
    const target = byId(kind === "downloader" ? "downloader-panel" : "update-panel");
    target?.scrollIntoView({ behavior: "smooth", block: "start" });
    PhaseUI.downloadMap?.queueRender();
  }, 30);
}

function renderSlaveInventory() {
  const isSEN = constellationCode() === "SEN";
  const files = asArray(PhaseUI.state?.slaves).filter(file => isSEN
    ? String(file?.type || "") === "Sentinel-1"
    : ["CSK", "CSG"].includes(String(file?.type || "")));
  byId("slaves-count").textContent = `${files.length} file${files.length === 1 ? "" : "s"}`;
  const body = byId("slaves-table-body");
  if (!files.length) {
    body.innerHTML = `<tr><td colspan="5" class="empty-cell">No ${isSEN ? "Sentinel-1 ZIP" : "COSMO-SkyMed HDF5"} files in slaves</td></tr>`;
    return;
  }
  body.innerHTML = files.map(file => `<tr>
    <td class="file-cell" title="${escapeAttribute(file.relativePath || file.name)}">${escapeHtml(file.name)}</td>
    <td>${escapeHtml(file.date || "—")}</td><td>${escapeHtml(file.type || "—")}</td>
    <td class="${/ready|processed/i.test(file.status || "") ? "status-ready" : ""}">${escapeHtml(file.status || "—")}</td>
    <td>${formatBytes(file.sizeBytes)}</td></tr>`).join("");
}

function renderDownloader() {
  const downloader = PhaseUI.state?.downloader || {};
  const results = asArray(downloader.results);
  const signature = results.map(item => item.sceneName).join("|");
  const resultsChanged = signature !== PhaseUI.asfResultSignature;
  if (resultsChanged) {
    PhaseUI.asfResultSignature = signature;
    PhaseUI.selectedAsfNames = new Set(results.filter(item => item.selected).map(item => String(item.sceneName)));
  }
  if (!PhaseUI.asfFiltersLoaded) {
    populateAsfFilters(downloader.filters || defaultAsfFilters());
    PhaseUI.asfFiltersLoaded = true;
  }

  const loggedIn = Boolean(downloader.loggedIn);
  byId("earthdata-signed-out").classList.toggle("hidden", loggedIn);
  byId("earthdata-signed-in").classList.toggle("hidden", !loggedIn);
  byId("earthdata-user").textContent = loggedIn ? `Signed in as ${downloader.username || "Earthdata user"}` : "Signed in";
  if (!loggedIn && downloader.username && !byId("earthdata-username").value) {
    byId("earthdata-username").value = downloader.username;
  }

  const sortedResults = sortAsfResults(results);
  renderAsfResults(sortedResults, Boolean(downloader.busy));
  renderDownloadMap(results, resultsChanged);
  byId("asf-recommended").textContent = downloader.recommended || "Run a search to load compatible acquisitions.";
  byId("asf-result-summary").textContent = `${results.length} product${results.length === 1 ? "" : "s"} · ${formatGigabytes(downloader.totalSizeGB)}`;
  byId("downloader-status").textContent = downloader.status || "Downloader ready";
  byId("asf-sort-key").value = PhaseUI.asfSortKey;
  byId("asf-sort-direction").textContent = sortDirectionLabel(
    PhaseUI.asfSortKey, PhaseUI.asfSortDirection
  );
  byId("search-asf").disabled = Boolean(downloader.busy);
  byId("earthdata-login").disabled = Boolean(downloader.busy);
}

function sortAsfResults(results) {
  const direction = PhaseUI.asfSortDirection === "desc" ? -1 : 1;
  const key = PhaseUI.asfSortKey;
  return [...results].sort((left, right) => {
    let comparison;
    if (key === "path") comparison = Number(left.pathNumber || 0) - Number(right.pathNumber || 0);
    else if (key === "frame") comparison = Number(left.frameNumber || 0) - Number(right.frameNumber || 0);
    else comparison = String(left.date || "").localeCompare(String(right.date || ""));
    if (!comparison) comparison = String(left.sceneName || "").localeCompare(String(right.sceneName || ""));
    return comparison * direction;
  });
}

function sortDirectionLabel(key, direction) {
  if (key === "date") return direction === "asc" ? "Oldest first" : "Newest first";
  return direction === "asc" ? "Lowest first" : "Highest first";
}

function renderDownloadMap(results, fitResults) {
  if (!PhaseUI.downloadMap) return;
  const selected = PhaseUI.selectedAsfNames;
  const footprints = results.map((item, index) => ({
    id: item.id || `asf-${index + 1}`,
    name: item.sceneName || `ASF product ${index + 1}`,
    source: "ASF search",
    selected: selected.has(String(item.sceneName)),
    coordinates: matrix(item.coordinates),
  }));
  const polygon = matrix(PhaseUI.state?.downloader?.polygon);
  PhaseUI.downloadMap.setData({
    coastlines: asArray(PhaseUI.state?.map?.coastlines),
    footprints,
    polygon,
  });
  updateDownloadAoiSummary(polygon);
  if (fitResults && footprints.length && PhaseUI.downloaderOpen) {
    window.setTimeout(() => PhaseUI.downloadMap.fitToFootprints(false), 20);
  }
}

function renderAsfResults(results, busy) {
  const body = byId("asf-results-body");
  if (!results.length) {
    body.innerHTML = `<tr><td colspan="8" class="empty-cell">No ASF results</td></tr>`;
  } else {
    body.innerHTML = results.map(item => {
      const name = String(item.sceneName || "");
      return `<tr><td><input class="asf-result-check" type="checkbox" data-scene="${escapeAttribute(name)}" ${PhaseUI.selectedAsfNames.has(name) ? "checked" : ""}></td>
        <td>${escapeHtml(item.date || "—")}</td><td>${escapeHtml(item.time || "—")}</td>
        <td>${escapeHtml(item.pathNumber)}</td><td>${escapeHtml(item.frameNumber)}</td>
        <td>${escapeHtml(item.direction || "—")}</td><td>${formatGigabytes(item.sizeGB)}</td>
        <td class="file-cell" title="${escapeAttribute(name)}">${escapeHtml(name)}</td></tr>`;
    }).join("");
    body.querySelectorAll(".asf-result-check").forEach(control => control.addEventListener("change", () => {
      if (control.checked) PhaseUI.selectedAsfNames.add(control.dataset.scene);
      else PhaseUI.selectedAsfNames.delete(control.dataset.scene);
      updateAsfSelectionControls(results, busy);
      renderDownloadMap(results, false);
    }));
  }
  updateAsfSelectionControls(results, busy);
}

function updateAsfSelectionControls(results, busy) {
  const selectedCount = results.filter(item => PhaseUI.selectedAsfNames.has(String(item.sceneName))).length;
  const selectAll = byId("asf-select-all");
  selectAll.checked = Boolean(results.length) && selectedCount === results.length;
  selectAll.indeterminate = selectedCount > 0 && selectedCount < results.length;
  selectAll.disabled = busy || !results.length;
  const button = byId("download-selected");
  button.disabled = busy || selectedCount === 0;
  button.textContent = selectedCount ? `Download selected (${selectedCount})` : "Download selected";
}

function selectAllAsf(checked) {
  const results = asArray(PhaseUI.state?.downloader?.results);
  PhaseUI.selectedAsfNames = new Set(checked ? results.map(item => String(item.sceneName)) : []);
  renderAsfResults(results, Boolean(PhaseUI.state?.downloader?.busy));
  renderDownloadMap(results, false);
}

function handleDownloadPolygon(polygon) {
  if (!polygon || polygon.length < 3) return;
  const closed = closePolygon(polygon);
  if (PhaseUI.state?.downloader) PhaseUI.state.downloader.polygon = closed;
  updateDownloadAoiSummary(closed);
  send("DownloadAoiChanged", { polygon: closed });
}

function useProcessingAoiForDownload() {
  const polygon = matrix(PhaseUI.state?.map?.polygon);
  if (polygon.length < 4) {
    toast("Draw the processing AOI first, then reuse it here.");
    return;
  }
  PhaseUI.downloadMap.setPolygon(polygon, true);
  handleDownloadPolygon(polygon);
}

function updateDownloadAoiSummary(polygon) {
  const closed = matrix(polygon);
  if (closed.length < 3) {
    byId("download-aoi-summary").textContent = "No download AOI selected";
    return;
  }
  const lons = closed.map(point => point[0]);
  const lats = closed.map(point => point[1]);
  byId("download-aoi-summary").textContent = `AOI ${formatCoordinate(Math.min(...lons))}, ${formatCoordinate(Math.min(...lats))} → ${formatCoordinate(Math.max(...lons))}, ${formatCoordinate(Math.max(...lats))}`;
}

function defaultAsfFilters() {
  return {
    dataset: "SENTINEL-1", processingLevel: ["SLC"], beamMode: ["IW"],
    polarization: [], flightDirection: [], subtype: [], startDate: "", endDate: "",
    pathStart: "", pathEnd: "", frameStart: "", frameEnd: "", groupID: "",
    samplingRate: "", samplingUnit: "Month",
  };
}

function populateAsfFilters(filters) {
  const values = { ...defaultAsfFilters(), ...(filters || {}) };
  setValue("asf-start-date", dateInputValue(values.startDate));
  setValue("asf-end-date", dateInputValue(values.endDate));
  setValue("asf-path-start", values.pathStart);
  setValue("asf-path-end", values.pathEnd);
  setValue("asf-frame-start", values.frameStart);
  setValue("asf-frame-end", values.frameEnd);
  setValue("asf-group-id", values.groupID);
  setValue("asf-sampling-rate", values.samplingRate);
  setValue("asf-sampling-unit", values.samplingUnit || "Month");
  setCheckedValues("filter-processing-level", values.processingLevel);
  setCheckedValues("filter-beam-mode", values.beamMode);
  setCheckedValues("filter-polarization", values.polarization);
  setCheckedValues("filter-flight-direction", values.flightDirection);
  setCheckedValues("filter-subtype", values.subtype);
}

function collectAsfFilters() {
  return {
    dataset: "SENTINEL-1",
    processingLevel: checkedValues("filter-processing-level"),
    beamMode: checkedValues("filter-beam-mode"),
    polarization: checkedValues("filter-polarization"),
    flightDirection: checkedValues("filter-flight-direction"),
    subtype: checkedValues("filter-subtype"),
    startDate: byId("asf-start-date").value,
    endDate: byId("asf-end-date").value,
    pathStart: byId("asf-path-start").value,
    pathEnd: byId("asf-path-end").value,
    frameStart: byId("asf-frame-start").value,
    frameEnd: byId("asf-frame-end").value,
    groupID: byId("asf-group-id").value,
    samplingRate: byId("asf-sampling-rate").value,
    samplingUnit: byId("asf-sampling-unit").value,
  };
}

function resetAsfFilters() {
  populateAsfFilters(defaultAsfFilters());
  toast("ASF filters reset to Sentinel-1 SLC / IW.");
}

function searchAsf() {
  collectVisibleForm();
  send("DownloadSearch", { ...payload(), filters: collectAsfFilters() });
}

function loginEarthdata() {
  const username = byId("earthdata-username").value.trim();
  const password = byId("earthdata-password").value;
  if (!username || !password) {
    toast("Enter both Earthdata username and password.");
    return;
  }
  send("DownloadLogin", { username, password });
  byId("earthdata-password").value = "";
}

function downloadSelectedAsf() {
  collectVisibleForm();
  send("DownloadSelected", { ...payload(), sceneNames: [...PhaseUI.selectedAsfNames] });
}

function renderUpdate() {
  const update = PhaseUI.state?.update || {};
  const context = update.context || {};
  const results = asArray(update.results);
  const signature = results.map(item => item.sceneName).join("|");
  if (signature !== PhaseUI.updateResultSignature) {
    PhaseUI.updateResultSignature = signature;
    PhaseUI.selectedUpdateNames = new Set(results.filter(item => item.selected !== false).map(item => String(item.sceneName)));
  }
  byId("update-context").textContent = context.message || "No local Sentinel-1 stack context is available.";
  if (!byId("update-end-date").value) byId("update-end-date").value = todayInputValue();
  byId("search-update").disabled = Boolean(update.busy) || !context.available;
  byId("refresh-update").disabled = Boolean(update.busy);
  byId("update-status").textContent = update.status || "Stack update ready";
  byId("update-sort-key").value = PhaseUI.updateSortKey;
  byId("update-sort-direction").textContent = sortDirectionLabel(
    PhaseUI.updateSortKey, PhaseUI.updateSortDirection
  );
  renderUpdateResults(sortUpdateResults(results), Boolean(update.busy));
}

function sortUpdateResults(results) {
  const direction = PhaseUI.updateSortDirection === "desc" ? -1 : 1;
  const key = PhaseUI.updateSortKey;
  return [...results].sort((left, right) => {
    let comparison;
    if (key === "path") comparison = Number(left.pathNumber || 0) - Number(right.pathNumber || 0);
    else if (key === "frame") comparison = Number(left.frameNumber || 0) - Number(right.frameNumber || 0);
    else comparison = String(left.date || left.startTime || "")
      .localeCompare(String(right.date || right.startTime || ""));
    return (comparison || String(left.sceneName || "").localeCompare(String(right.sceneName || ""))) * direction;
  });
}

function renderUpdateResults(results, busy) {
  const body = byId("update-results-body");
  if (!results.length) {
    body.innerHTML = `<tr><td colspan="7" class="empty-cell">No update search results</td></tr>`;
  } else {
    body.innerHTML = results.map(item => {
      const name = String(item.sceneName || "");
      return `<tr><td><input class="update-result-check" type="checkbox" data-scene="${escapeAttribute(name)}" ${PhaseUI.selectedUpdateNames.has(name) ? "checked" : ""}></td>
        <td>${escapeHtml(item.date || item.startTime || "—")}</td>
        <td>${escapeHtml(item.pathNumber || "—")}</td><td>${escapeHtml(item.frameNumber || "—")}</td>
        <td>${escapeHtml(item.platform || "Sentinel-1")}</td>
        <td>${escapeHtml(item.polarization || "—")}</td><td class="file-cell" title="${escapeAttribute(name)}">${escapeHtml(name)}</td></tr>`;
    }).join("");
    body.querySelectorAll(".update-result-check").forEach(control => control.addEventListener("change", () => {
      if (control.checked) PhaseUI.selectedUpdateNames.add(control.dataset.scene);
      else PhaseUI.selectedUpdateNames.delete(control.dataset.scene);
      updateUpdateSelectionControls(results, busy);
    }));
  }
  updateUpdateSelectionControls(results, busy);
}

function updateUpdateSelectionControls(results, busy) {
  const selectedCount = results.filter(item => PhaseUI.selectedUpdateNames.has(String(item.sceneName))).length;
  const selectAll = byId("update-select-all");
  selectAll.checked = Boolean(results.length) && selectedCount === results.length;
  selectAll.indeterminate = selectedCount > 0 && selectedCount < results.length;
  selectAll.disabled = busy || !results.length;
  const button = byId("download-update");
  button.disabled = busy || selectedCount === 0;
  button.textContent = selectedCount ? `Download selected (${selectedCount})` : "Download selected";
}

function selectAllUpdate(checked) {
  const results = asArray(PhaseUI.state?.update?.results);
  PhaseUI.selectedUpdateNames = new Set(checked ? results.map(item => String(item.sceneName)) : []);
  renderUpdateResults(results, Boolean(PhaseUI.state?.update?.busy));
}

function searchUpdate() {
  send("SearchUpdate", { endDate: compactDate(byId("update-end-date").value) });
}

function downloadSelectedUpdate() {
  send("DownloadUpdate", { sceneNames: [...PhaseUI.selectedUpdateNames] });
}

function renderTransfer() {
  const transfer = PhaseUI.state?.transfer || {};
  renderTransferCard("download", transfer, transfer.kind === "initial");
  renderTransferCard("update", transfer, transfer.kind === "update");
}

function renderTransferCard(prefix, transfer, matches) {
  const card = byId(`${prefix}-transfer-card`);
  if (!card) return;
  card.classList.toggle("hidden", !matches);
  if (!matches) return;
  const percentage = Math.max(0, Math.min(100, Number(transfer.percentage) || 0));
  const currentIndex = Number(transfer.currentIndex) || 0;
  const totalFiles = Number(transfer.totalFiles) || 0;
  const currentBytes = Number(transfer.currentBytes) || 0;
  const currentTotal = Number(transfer.currentTotalBytes) || 0;
  byId(`${prefix}-transfer-percent`).textContent = `${percentage.toFixed(1)}%`;
  byId(`${prefix}-transfer-progress`).style.width = `${percentage}%`;
  byId(`${prefix}-transfer-count`).textContent = `Image ${currentIndex} of ${totalFiles}`;
  byId(`${prefix}-transfer-bytes`).textContent = currentTotal
    ? `${formatBytes(currentBytes)} of ${formatBytes(currentTotal)}`
    : `${formatBytes(currentBytes)} downloaded`;
  byId(`${prefix}-transfer-file`).textContent = transfer.currentFile || transfer.message || "Waiting to start";
  const messageNode = prefix === "download" ? byId("downloader-status") : byId("update-transfer-message");
  messageNode.textContent = transfer.message || "Download ready";
  const stopButton = byId(`stop-${prefix}-transfer`);
  stopButton.disabled = !transfer.canStop;
  card.classList.toggle("complete", transfer.phase === "completed");
  card.classList.toggle("failed", transfer.phase === "failed");
  card.classList.toggle("stopped", transfer.phase === "stopped");
}

function appendLog(entry) {
  if (!entry) return;
  PhaseUI.logs.push(entry);
  const consoleNode = byId("console");
  if (consoleNode) { consoleNode.appendChild(logElement(entry)); consoleNode.scrollTop = consoleNode.scrollHeight; updateLogCount(); }
}

function renderLogs() {
  const consoleNode = byId("console");
  consoleNode.innerHTML = "";
  PhaseUI.logs.forEach(entry => consoleNode.appendChild(logElement(entry)));
  consoleNode.scrollTop = consoleNode.scrollHeight;
  updateLogCount();
}

function logElement(entry) {
  const line = document.createElement("div");
  line.className = "log-line";
  const message = typeof entry === "string" ? entry : (entry.message || "");
  const time = typeof entry === "string" ? "" : (entry.time || "");
  let kind = "";
  if (/error|failed|missing/i.test(message)) kind = "error";
  else if (/finished|completed|success|ready/i.test(message)) kind = "success";
  else if (/---+.*step/i.test(message)) kind = "step";
  line.innerHTML = `<span class="log-time">${escapeHtml(time)}</span><span class="log-message ${kind}">${escapeHtml(message)}</span>`;
  return line;
}

function updateLogCount() { byId("log-count").textContent = `${PhaseUI.logs.length} message${PhaseUI.logs.length === 1 ? "" : "s"}`; }
function save() { collectVisibleForm(); send("Save", payload()); }
function start() { collectVisibleForm(); PhaseUI.active = "run"; renderNavigation(); renderPage(); send("Start", payload()); }
function payload() { return { config: { ...PhaseUI.config } }; }
function send(name, data) {
  if (!PhaseUI.bridge) { toast("MATLAB bridge is not connected."); return; }
  PhaseUI.bridge.sendEventToMATLAB(name, data || {});
}
function moveSection(delta) {
  collectVisibleForm();
  const groups = PhaseUI.schema.groups;
  const index = groups.findIndex(group => group.id === PhaseUI.active);
  const next = groups[index + delta];
  if (next) switchPage(next.id);
}
function constellationCode() { return /^COSMO/i.test(PhaseUI.config.constellation || "") ? "CSK" : "SEN"; }
function bboxPolygonFromConfig() {
  const minLon = Number(PhaseUI.config.lon_min), maxLon = Number(PhaseUI.config.lon_max);
  const minLat = Number(PhaseUI.config.lat_min), maxLat = Number(PhaseUI.config.lat_max);
  return [[minLon,minLat],[maxLon,minLat],[maxLon,maxLat],[minLon,maxLat],[minLon,minLat]];
}
function matrix(value) {
  if (!value) return [];
  if (Array.isArray(value) && value.length && Array.isArray(value[0])) return value.map(row => [Number(row[0]), Number(row[1])]);
  return [];
}
function closePolygon(points) {
  const result = points.map(point => [Number(point[0]), Number(point[1])]);
  if (!result.length) return result;
  const first = result[0], last = result[result.length - 1];
  if (first[0] !== last[0] || first[1] !== last[1]) result.push([...first]);
  return result;
}
function formatCoordinate(value) { return Number(value).toFixed(6).replace(/0+$/, "").replace(/\.$/, ""); }
function findGroup(id) { return PhaseUI.schema.groups.find(group => group.id === id); }
function byId(id) { return document.getElementById(id); }
function listen(id, eventName, callback) {
  const node = byId(id);
  if (node) node.addEventListener(eventName, callback);
}
function asArray(value) { if (!value) return []; return Array.isArray(value) ? value : [value]; }
function textValue(value) { return value === null || value === undefined ? "" : String(value); }
function escapeHtml(value) { const node = document.createElement("div"); node.textContent = textValue(value); return node.innerHTML; }
function escapeAttribute(value) { return escapeHtml(value).replace(/`/g, "&#96;"); }
function checkedValues(containerId) {
  return [...byId(containerId).querySelectorAll('input[type="checkbox"]:checked')].map(input => input.value);
}
function setCheckedValues(containerId, values) {
  const selected = new Set(asArray(values).map(String));
  byId(containerId).querySelectorAll('input[type="checkbox"]').forEach(input => { input.checked = selected.has(input.value); });
}
function setValue(id, value) { byId(id).value = value === null || value === undefined ? "" : String(value); }
function compactDate(value) { return String(value || "").replace(/[^0-9]/g, ""); }
function dateInputValue(value) {
  const raw = compactDate(value);
  return raw.length === 8 ? `${raw.slice(0,4)}-${raw.slice(4,6)}-${raw.slice(6,8)}` : "";
}
function todayInputValue() {
  const now = new Date();
  const local = new Date(now.getTime() - now.getTimezoneOffset() * 60000);
  return local.toISOString().slice(0, 10);
}
function formatBytes(value) {
  const bytes = Number(value) || 0;
  if (bytes >= 1024 ** 3) return `${(bytes / 1024 ** 3).toFixed(2)} GB`;
  if (bytes >= 1024 ** 2) return `${(bytes / 1024 ** 2).toFixed(1)} MB`;
  if (bytes >= 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  return `${bytes} B`;
}
function formatGigabytes(value) { return `${(Number(value) || 0).toFixed(2)} GB`; }
function toast(message) {
  const node = document.createElement("div"); node.className = "toast"; node.textContent = message;
  byId("toast-stack").appendChild(node); setTimeout(() => node.remove(), 3600);
}

window.addEventListener("error", event => {
  const banner = byId("error-banner");
  if (!banner) return;
  banner.textContent = `Interface error: ${event.message}`; banner.classList.remove("hidden");
});
