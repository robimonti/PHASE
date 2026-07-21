"use strict";

const PhaseUI = {
  bridge: null,
  state: null,
  config: {},
  schema: { groups: [], items: [] },
  active: "project",
  localDirty: false,
  logs: [],
  detected: new Set(),
};

function setup(htmlComponent) {
  PhaseUI.bridge = htmlComponent;
  htmlComponent.addEventListener("DataChanged", () => receiveState(htmlComponent.Data));
  htmlComponent.addEventListener("PhaseLog", event => appendLog(event.Data));
  document.addEventListener("DOMContentLoaded", wireStaticControls, { once: true });
  if (document.readyState !== "loading") wireStaticControls();
  htmlComponent.sendEventToMATLAB("Ready", {});
  if (htmlComponent.Data) receiveState(htmlComponent.Data);
}

function wireStaticControls() {
  if (wireStaticControls.done) return;
  wireStaticControls.done = true;
  byId("load-button").addEventListener("click", () => send("Load", {}));
  byId("save-button").addEventListener("click", save);
  byId("banner-save").addEventListener("click", save);
  byId("start-button").addEventListener("click", start);
  byId("workdir").addEventListener("click", () => send("OpenWorkDir", {}));
  byId("open-error-log").addEventListener("click", () => send("OpenErrorLog", {}));
  byId("open-ts-picker").addEventListener("click", () => send("OpenTsPicker", payload()));
  byId("clear-log").addEventListener("click", () => { PhaseUI.logs = []; renderLogs(); });
  byId("advanced-toggle").addEventListener("change", event => {
    document.body.classList.toggle("show-advanced", event.target.checked);
  });
  byId("previous-section").addEventListener("click", () => moveSection(-1));
  byId("next-section").addEventListener("click", () => moveSection(1));
}

function receiveState(state) {
  if (!state || state.kind !== "state") return;
  PhaseUI.state = state;
  PhaseUI.schema = state.schema || { groups: [], items: [] };
  PhaseUI.config = { ...(state.config || {}) };
  PhaseUI.localDirty = Boolean(state.dirty);
  PhaseUI.detected = new Set(asArray(state.detectedFields));
  PhaseUI.logs = asArray(state.logs);
  if (!findGroup(PhaseUI.active)) PhaseUI.active = "project";
  renderAll();
}

function renderAll() {
  renderNavigation();
  renderPage();
  renderStatus();
  renderLogs();
  renderSummary();
  byId("workdir").textContent = PhaseUI.state?.workDir || "No processing folder";
  byId("version").textContent = `PHASE StaMPS ${PhaseUI.state?.version || "beta"}`;
  byId("blocking-overlay").classList.toggle("hidden", !PhaseUI.state?.running);
}

function renderNavigation() {
  const nav = byId("navigation");
  nav.innerHTML = `<div class="nav-label">Configuration</div>`;
  PhaseUI.schema.groups.forEach((group, index) => {
    const button = document.createElement("button");
    button.type = "button";
    button.className = `nav-item${PhaseUI.active === group.id ? " active" : ""}`;
    button.innerHTML = `<span class="nav-index">${index + 1}</span><span class="nav-text">${escapeHtml(group.title)}</span>${PhaseUI.localDirty ? '<i class="nav-dirty"></i>' : ""}`;
    button.addEventListener("click", () => switchPage(group.id));
    nav.appendChild(button);
  });
  nav.insertAdjacentHTML("beforeend", `<div class="nav-label">Tools</div>`);
  nav.appendChild(toolNav("run", "⌁", "Run monitor"));
  nav.appendChild(toolNav("ts", "⌖", "TS Points"));
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
  ["form-page", "run-page", "ts-page"].forEach(id => byId(id).classList.remove("active"));
  if (PhaseUI.active === "run") { byId("run-page").classList.add("active"); renderSummary(); return; }
  if (PhaseUI.active === "ts") { byId("ts-page").classList.add("active"); return; }
  byId("form-page").classList.add("active");
  const group = findGroup(PhaseUI.active) || PhaseUI.schema.groups[0];
  if (!group) return;
  byId("section-title").textContent = group.title;
  byId("section-subtitle").textContent = group.subtitle || "";
  const items = PhaseUI.schema.items.filter(item => item.group === group.id);
  byId("section-meta").textContent = `${items.length} controls`;
  const grid = byId("form-grid");
  grid.innerHTML = "";
  items.forEach(item => grid.appendChild(createField(item)));
  const index = PhaseUI.schema.groups.findIndex(candidate => candidate.id === group.id);
  byId("previous-section").disabled = index <= 0;
  byId("next-section").disabled = index >= PhaseUI.schema.groups.length - 1;
}

function createField(item) {
  const wrapper = document.createElement("div");
  wrapper.className = `field${item.advanced ? " advanced" : ""}${item.type === "path" ? " full" : ""}${PhaseUI.detected.has(item.id) ? " detected" : ""}`;
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
      element.value = String(option); element.textContent = String(option); control.appendChild(element);
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
    control.type = item.type === "date" ? "date" : "text";
    control.value = textValue(PhaseUI.config[item.id]);
    control.spellcheck = false;
  }
  control.id = `field-${item.id}`;
  control.dataset.configField = item.id;
  control.addEventListener("input", markDirty);
  control.addEventListener("change", markDirty);

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
}

function markDirty() {
  collectVisibleForm();
  PhaseUI.localDirty = true;
  renderStatus();
  renderNavigation();
  renderSummary();
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
  const first = Number(PhaseUI.config.stamps_first_step || 1);
  const last = Number(PhaseUI.config.stamps_last_step || 7);
  const width = Math.max(0, Math.min(100, ((last - first + 1) / 8) * 100));
  byId("range-progress").style.width = `${width}%`;
  byId("range-label").textContent = `StaMPS steps ${first} → ${last}`;
}

function renderSummary() {
  byId("summary-first").textContent = textValue(PhaseUI.config.stamps_first_step) || "—";
  byId("summary-last").textContent = textValue(PhaseUI.config.stamps_last_step) || "—";
  byId("summary-tropo").textContent = PhaseUI.config.train_enabled && PhaseUI.config.subtr_tropo === "y" ? `TRAIN · ${PhaseUI.config.tropo_method}` : "Disabled";
  byId("summary-output").textContent = PhaseUI.config.ph_output || "—";
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
  else if (/finished|completed|success/i.test(message)) kind = "success";
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

function findGroup(id) { return PhaseUI.schema.groups.find(group => group.id === id); }
function byId(id) { return document.getElementById(id); }
function asArray(value) { if (!value) return []; return Array.isArray(value) ? value : [value]; }
function textValue(value) { return value === null || value === undefined ? "" : String(value); }
function escapeHtml(value) { const node = document.createElement("div"); node.textContent = textValue(value); return node.innerHTML; }
function toast(message) {
  const node = document.createElement("div"); node.className = "toast"; node.textContent = message;
  byId("toast-stack").appendChild(node); setTimeout(() => node.remove(), 3600);
}

window.addEventListener("error", event => {
  const banner = byId("error-banner");
  if (!banner) return;
  banner.textContent = `Interface error: ${event.message}`; banner.classList.remove("hidden");
});
