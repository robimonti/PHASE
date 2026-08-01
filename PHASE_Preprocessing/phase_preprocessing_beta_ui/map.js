"use strict";

/*
 * Dependency-free vector map used inside MATLAB uihtml.
 * Coordinates crossing the MATLAB/JavaScript boundary are always [lon, lat].
 */
class PhaseMap {
  constructor(container, options = {}) {
    this.container = container;
    this.svg = container.querySelector("svg");
    this.tooltip = container.querySelector(".map-tooltip");
    this.coordinate = container.querySelector(".map-coordinate");
    this.onPolygonChanged = options.onPolygonChanged || (() => {});
    this.onDrawingChanged = options.onDrawingChanged || (() => {});
    this.onTilesRequested = options.onTilesRequested || (() => {});

    this.coastlines = [];
    this.footprints = [];
    this.polygon = [];
    this.draft = [];
    this.cursorPoint = null;
    this.drawing = false;
    this.initialised = false;
    this.center = { x: 0.535, y: 0.375 };
    this.zoom = 4;
    this.pan = null;
    this.vertexDrag = null;
    this.resizeFrame = 0;
    this.tileRequestTimer = 0;
    this.pendingTileRequests = new Map();
    this.requestedTiles = new Set();

    this.imageryLayer = document.createElement("div");
    this.imageryLayer.className = "map-tiles map-imagery-tiles";
    this.labelLayer = document.createElement("div");
    this.labelLayer.className = "map-tiles map-label-tiles";
    this.container.insertBefore(this.imageryLayer, this.svg);
    this.container.insertBefore(this.labelLayer, this.svg);

    this.layers = {};
    ["grid", "coast", "footprints", "aoi", "vertices"].forEach(name => {
      const group = svgNode("g", { class: `map-layer map-${name}-layer` });
      this.svg.appendChild(group);
      this.layers[name] = group;
    });

    this.bindInteractions();
    if (typeof ResizeObserver !== "undefined") {
      this.resizeObserver = new ResizeObserver(() => this.queueRender());
      this.resizeObserver.observe(this.container);
    }
    window.addEventListener("resize", () => this.queueRender());
  }

  setData(data = {}) {
    this.coastlines = asCollection(data.coastlines).map(item => normaliseMatrix(item?.coordinates));
    this.footprints = asCollection(data.footprints).map((item, index) => ({
      id: String(item?.id || `footprint-${index + 1}`),
      name: String(item?.name || `Footprint ${index + 1}`),
      source: String(item?.source || "Satellite product"),
      selected: Boolean(item?.selected),
      coordinates: normaliseMatrix(item?.coordinates),
    })).filter(item => item.coordinates.length >= 3);

    if (!this.drawing && !this.vertexDrag) this.polygon = closePolygon(normaliseMatrix(data.polygon));
    if (!this.initialised && this.container.clientWidth > 40 && this.container.clientHeight > 40) {
      this.initialised = true;
      if (this.polygon.length) this.fitBounds(this.polygon, false);
      else if (this.footprints.length) this.fitToFootprints(false);
    }
    this.render();
  }

  setPolygon(polygon, fit = false) {
    this.polygon = closePolygon(normaliseMatrix(polygon));
    this.draft = [];
    this.cursorPoint = null;
    if (fit && this.polygon.length) this.fitBounds(this.polygon, false);
    this.render();
  }

  startDrawing() {
    this.drawing = true;
    this.draft = [];
    this.cursorPoint = null;
    this.hideTooltip();
    this.container.classList.add("drawing");
    this.container.focus();
    this.onDrawingChanged(true);
    this.updateCoordinate(null, "Click to add the first polygon vertex · Enter to finish · Esc to cancel");
    this.render();
  }

  finishDrawing() {
    const polygon = closePolygon(removeRepeatedPoints(this.draft));
    if (polygon.length < 4) {
      this.updateCoordinate(null, `Add ${Math.max(0, 3 - this.draft.length)} more point${this.draft.length === 2 ? "" : "s"} to finish`);
      return;
    }
    this.drawing = false;
    this.draft = [];
    this.cursorPoint = null;
    this.polygon = polygon;
    this.container.classList.remove("drawing");
    this.onDrawingChanged(false);
    this.updateCoordinate(null, "Polygon saved · drag a vertex to refine it");
    this.render();
    this.onPolygonChanged(cloneMatrix(this.polygon));
  }

  cancelDrawing() {
    if (!this.drawing) return;
    this.drawing = false;
    this.draft = [];
    this.cursorPoint = null;
    this.container.classList.remove("drawing");
    this.onDrawingChanged(false);
    this.updateCoordinate(null, "Drawing cancelled");
    this.render();
  }

  fitToFootprints(animate = true) {
    const selected = this.footprints.filter(item => item.selected);
    const source = selected.length ? selected : this.footprints;
    const points = source.flatMap(item => item.coordinates);
    if (points.length) this.fitBounds(points, animate);
    else if (this.polygon.length) this.fitBounds(this.polygon, animate);
  }

  fitBounds(points, animate = true) {
    const projected = normaliseMatrix(points).map(point => project(point[0], point[1]));
    if (!projected.length) return;
    const width = Math.max(320, this.container.clientWidth || 800);
    const height = Math.max(260, this.container.clientHeight || 480);
    const xs = unwrapWorldXs(projected.map(point => point.x));
    const ys = projected.map(point => point.y);
    const minX = Math.min(...xs), maxX = Math.max(...xs);
    const minY = Math.min(...ys), maxY = Math.max(...ys);
    const spanX = Math.max(maxX - minX, 0.00015);
    const spanY = Math.max(maxY - minY, 0.00015);
    const padding = Math.min(82, Math.max(38, Math.min(width, height) * 0.12));
    const scale = Math.min((width - padding * 2) / spanX, (height - padding * 2) / spanY);
    this.center = { x: wrap01((minX + maxX) / 2), y: clamp((minY + maxY) / 2, 0.002, 0.998) };
    this.zoom = clamp(Math.log2(scale / 256), 0.45, 13);
    if (animate) {
      this.container.classList.add("map-moving");
      window.setTimeout(() => this.container.classList.remove("map-moving"), 220);
    }
    this.render();
  }

  bindInteractions() {
    this.container.addEventListener("pointerdown", event => this.pointerDown(event));
    this.container.addEventListener("pointermove", event => this.pointerMove(event));
    this.container.addEventListener("pointerup", event => this.pointerUp(event));
    this.container.addEventListener("pointercancel", event => this.pointerUp(event));
    this.container.addEventListener("pointerleave", () => {
      if (!this.pan && !this.vertexDrag) this.updateCoordinate(null);
    });
    this.container.addEventListener("click", event => this.mapClick(event));
    this.container.addEventListener("dblclick", event => {
      if (!this.drawing) return;
      event.preventDefault();
      event.stopPropagation();
      this.finishDrawing();
    });
    this.container.addEventListener("wheel", event => this.wheel(event), { passive: false });
    this.container.addEventListener("keydown", event => {
      if (event.key === "Enter" && this.drawing) { event.preventDefault(); this.finishDrawing(); }
      if (event.key === "Escape" && this.drawing) { event.preventDefault(); this.cancelDrawing(); }
    });
  }

  pointerDown(event) {
    if (event.button !== 0) return;
    const vertex = event.target.closest?.(".map-vertex");
    if (vertex && !this.drawing) {
      event.preventDefault();
      event.stopPropagation();
      this.vertexDrag = { index: Number(vertex.dataset.index), pointerId: event.pointerId };
      this.container.setPointerCapture?.(event.pointerId);
      this.container.classList.add("vertex-dragging");
      return;
    }
    if (this.drawing) return;
    this.pan = {
      x: event.clientX,
      y: event.clientY,
      centerX: this.center.x,
      centerY: this.center.y,
      pointerId: event.pointerId,
      moved: false,
    };
    this.container.setPointerCapture?.(event.pointerId);
    this.container.classList.add("panning");
    this.hideTooltip();
  }

  pointerMove(event) {
    const geo = this.geoAtEvent(event);
    if (this.vertexDrag) {
      const point = [roundCoordinate(geo.lon), roundCoordinate(geo.lat)];
      const unique = openPolygon(this.polygon);
      if (this.vertexDrag.index >= 0 && this.vertexDrag.index < unique.length) {
        unique[this.vertexDrag.index] = point;
        this.polygon = closePolygon(unique);
        this.render();
      }
      this.updateCoordinate(geo, "Dragging AOI vertex");
      return;
    }
    if (this.pan) {
      const scale = this.worldScale();
      const dx = event.clientX - this.pan.x;
      const dy = event.clientY - this.pan.y;
      this.pan.moved = this.pan.moved || Math.abs(dx) + Math.abs(dy) > 3;
      this.center.x = wrap01(this.pan.centerX - dx / scale);
      this.center.y = clamp(this.pan.centerY - dy / scale, 0.002, 0.998);
      this.render();
      this.updateCoordinate(geo);
      return;
    }
    if (this.drawing) {
      this.cursorPoint = [roundCoordinate(geo.lon), roundCoordinate(geo.lat)];
      this.renderAoi();
      this.updateCoordinate(geo, `${this.draft.length} ${this.draft.length === 1 ? "vertex" : "vertices"} · click to add`);
    } else {
      this.updateCoordinate(geo);
    }
  }

  pointerUp(event) {
    if (this.vertexDrag) {
      this.container.releasePointerCapture?.(this.vertexDrag.pointerId);
      this.vertexDrag = null;
      this.container.classList.remove("vertex-dragging");
      this.onPolygonChanged(cloneMatrix(this.polygon));
      return;
    }
    if (!this.pan) return;
    this.container.releasePointerCapture?.(this.pan.pointerId);
    this.suppressClick = this.pan.moved;
    this.pan = null;
    this.container.classList.remove("panning");
  }

  mapClick(event) {
    if (this.suppressClick) { this.suppressClick = false; return; }
    if (!this.drawing || event.target.closest?.(".map-vertex")) return;
    const geo = this.geoAtEvent(event);
    const point = [roundCoordinate(geo.lon), roundCoordinate(geo.lat)];
    const previous = this.draft[this.draft.length - 1];
    if (!previous || Math.abs(previous[0] - point[0]) + Math.abs(previous[1] - point[1]) > 1e-7) {
      this.draft.push(point);
    }
    this.render();
    this.updateCoordinate(geo, `${this.draft.length} ${this.draft.length === 1 ? "vertex" : "vertices"} · double-click or Enter to finish`);
  }

  wheel(event) {
    event.preventDefault();
    const rect = this.container.getBoundingClientRect();
    const px = event.clientX - rect.left;
    const py = event.clientY - rect.top;
    const before = this.worldAtScreen(px, py);
    const delta = clamp(-event.deltaY * 0.0018, -0.9, 0.9);
    const nextZoom = clamp(this.zoom + delta, 0.45, 13);
    if (nextZoom === this.zoom) return;
    this.zoom = nextZoom;
    const scale = this.worldScale();
    this.center.x = wrap01(before.x - (px - rect.width / 2) / scale);
    this.center.y = clamp(before.y - (py - rect.height / 2) / scale, 0.002, 0.998);
    this.render();
  }

  render() {
    const width = Math.max(1, this.container.clientWidth);
    const height = Math.max(1, this.container.clientHeight);
    this.renderBasemap(width, height);
    this.svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
    this.renderGrid();
    this.renderCoastlines();
    this.renderFootprints();
    this.renderAoi();
  }

  renderBasemap(width, height) {
    if (width < 40 || height < 40) return;
    const tileZoom = clamp(Math.floor(this.zoom), 0, 18);
    const tileCount = 2 ** tileZoom;
    const tileSize = 256 * (2 ** (this.zoom - tileZoom));
    const centerTileX = this.center.x * tileCount;
    const centerTileY = this.center.y * tileCount;
    const minimumX = Math.floor(centerTileX - width / (2 * tileSize)) - 1;
    const maximumX = Math.ceil(centerTileX + width / (2 * tileSize)) + 1;
    const minimumY = Math.max(0, Math.floor(centerTileY - height / (2 * tileSize)) - 1);
    const maximumY = Math.min(tileCount - 1, Math.ceil(centerTileY + height / (2 * tileSize)) + 1);
    const layout = { tileZoom, tileCount, tileSize, centerTileX, centerTileY, minimumX, maximumX, minimumY, maximumY, width, height };
    this.renderTileLayer(this.imageryLayer, layout, "imagery");
    this.renderTileLayer(this.labelLayer, layout, "labels");
  }

  renderTileLayer(layer, layout, layerName) {
    if (!layer._tileNodes) layer._tileNodes = new Map();
    const visible = new Set();
    for (let y = layout.minimumY; y <= layout.maximumY; y += 1) {
      for (let rawX = layout.minimumX; rawX <= layout.maximumX; rawX += 1) {
        const x = ((rawX % layout.tileCount) + layout.tileCount) % layout.tileCount;
        const key = `${layout.tileZoom}/${rawX}/${y}`;
        visible.add(key);
        let tile = layer._tileNodes.get(key);
        if (!tile) {
          tile = document.createElement("img");
          tile.alt = "";
          tile.draggable = false;
          tile.decoding = "async";
          tile.dataset.cacheKey = `${layerName}/${layout.tileZoom}/${x}/${y}`;
          tile.dataset.localSource = localTileUrl(layerName, layout.tileZoom, x, y);
          tile.src = tile.dataset.localSource;
          tile.addEventListener("load", () => {
            tile.classList.add("loaded");
            tile.classList.remove("unavailable");
          });
          tile.addEventListener("error", () => {
            tile.classList.remove("loaded");
            tile.classList.add("unavailable");
            this.queueTileRequest({layer: layerName, z: layout.tileZoom, x, y, key: tile.dataset.cacheKey});
          });
          layer._tileNodes.set(key, tile);
          layer.appendChild(tile);
        }
        tile.style.width = `${layout.tileSize + 0.6}px`;
        tile.style.height = `${layout.tileSize + 0.6}px`;
        tile.style.left = `${layout.width / 2 + (rawX - layout.centerTileX) * layout.tileSize}px`;
        tile.style.top = `${layout.height / 2 + (y - layout.centerTileY) * layout.tileSize}px`;
      }
    }
    layer._tileNodes.forEach((tile, key) => {
      if (!visible.has(key)) {
        tile.remove();
        layer._tileNodes.delete(key);
      }
    });
  }

  queueTileRequest(request) {
    if (this.requestedTiles.has(request.key)) return;
    this.pendingTileRequests.set(request.key, request);
    window.clearTimeout(this.tileRequestTimer);
    this.tileRequestTimer = window.setTimeout(() => {
      const requests = [...this.pendingTileRequests.values()];
      this.pendingTileRequests.clear();
      requests.forEach(item => this.requestedTiles.add(item.key));
      if (requests.length) this.onTilesRequested(requests);
    }, 90);
  }

  retryTiles(keys) {
    const ready = new Set(asCollection(keys).map(String));
    if (!ready.size) return;
    ready.forEach(key => this.requestedTiles.delete(key));
    [this.imageryLayer, this.labelLayer].forEach(layer => {
      layer._tileNodes?.forEach(tile => {
        if (!ready.has(tile.dataset.cacheKey)) return;
        tile.classList.remove("loaded", "unavailable");
        tile.src = `${tile.dataset.localSource}?v=${Date.now()}`;
      });
    });
  }

  renderGrid() {
    clearNode(this.layers.grid);
    const interval = gridInterval(this.zoom);
    const scale = this.worldScale();
    const centreGeo = unproject(this.center.x, this.center.y);
    const longitudeSpan = Math.min(360, this.container.clientWidth / scale * 360);
    const minimumLongitude = centreGeo.lon - longitudeSpan / 2 - interval;
    const maximumLongitude = centreGeo.lon + longitudeSpan / 2 + interval;
    const topLatitude = unproject(this.center.x, this.center.y - this.container.clientHeight / scale / 2).lat;
    const bottomLatitude = unproject(this.center.x, this.center.y + this.container.clientHeight / scale / 2).lat;
    const minimumLatitude = clamp(bottomLatitude - interval, -85, 85);
    const maximumLatitude = clamp(topLatitude + interval, -85, 85);

    for (let rawLon = Math.ceil(minimumLongitude / interval) * interval; rawLon <= maximumLongitude; rawLon += interval) {
      const lon = wrapLongitude(rawLon);
      this.layers.grid.appendChild(this.pathNode([[lon, minimumLatitude], [lon, maximumLatitude]], "map-graticule"));
    }
    const latInterval = interval >= 30 ? 30 : interval >= 10 ? 10 : interval >= 2 ? 2 : interval;
    for (let lat = Math.ceil(minimumLatitude / latInterval) * latInterval; lat <= maximumLatitude; lat += latInterval) {
      const points = [];
      const samples = Math.max(2, Math.ceil(longitudeSpan / Math.max(interval, 1)));
      for (let index = 0; index <= samples; index += 1) {
        const rawLon = minimumLongitude + (maximumLongitude - minimumLongitude) * index / samples;
        points.push([wrapLongitude(rawLon), lat]);
      }
      this.layers.grid.appendChild(this.pathNode(points, lat === 0 ? "map-graticule equator" : "map-graticule"));
    }
  }

  renderCoastlines() {
    clearNode(this.layers.coast);
    this.coastlines.forEach(line => {
      if (line.length > 1) this.layers.coast.appendChild(this.pathNode(line, "map-coastline"));
    });
  }

  renderFootprints() {
    clearNode(this.layers.footprints);
    this.footprints.forEach(item => {
      const path = this.pathNode(closePolygon(item.coordinates), `map-footprint${item.selected ? " selected" : ""}`, true);
      path.setAttribute("tabindex", "0");
      path.setAttribute("aria-label", `${item.name}, ${item.source}`);
      path.addEventListener("pointerenter", event => this.showTooltip(item, event));
      path.addEventListener("pointermove", event => this.positionTooltip(event));
      path.addEventListener("pointerleave", () => this.hideTooltip());
      path.addEventListener("focus", event => this.showTooltip(item, event));
      path.addEventListener("blur", () => this.hideTooltip());
      this.layers.footprints.appendChild(path);
    });
  }

  renderAoi() {
    clearNode(this.layers.aoi);
    clearNode(this.layers.vertices);
    if (this.drawing) {
      const preview = [...this.draft];
      if (this.cursorPoint) preview.push(this.cursorPoint);
      if (preview.length) this.layers.aoi.appendChild(this.pathNode(preview, "map-aoi draft", false));
      this.draft.forEach((point, index) => this.layers.vertices.appendChild(this.vertexNode(point, index, true)));
      return;
    }
    if (this.polygon.length >= 4) {
      this.layers.aoi.appendChild(this.pathNode(this.polygon, "map-aoi", true));
      openPolygon(this.polygon).forEach((point, index) => {
        this.layers.vertices.appendChild(this.vertexNode(point, index, false));
      });
    }
  }

  pathNode(points, className, close = false) {
    const commands = [];
    let previous = null;
    normaliseMatrix(points).forEach(point => {
      const screen = this.screenAtGeo(point[0], point[1]);
      const farJump = previous && Math.abs(screen.x - previous.x) > Math.max(this.container.clientWidth * 0.72, 420);
      commands.push(`${!previous || farJump ? "M" : "L"}${screen.x.toFixed(2)},${screen.y.toFixed(2)}`);
      previous = screen;
    });
    if (close && commands.length) commands.push("Z");
    return svgNode("path", { d: commands.join(" "), class: className });
  }

  vertexNode(point, index, draft) {
    const screen = this.screenAtGeo(point[0], point[1]);
    return svgNode("circle", {
      cx: screen.x.toFixed(2), cy: screen.y.toFixed(2), r: draft ? 4.3 : 5,
      class: `map-vertex${draft ? " draft" : ""}`, "data-index": String(index),
    });
  }

  showTooltip(item, event) {
    if (this.drawing || this.pan) return;
    this.tooltip.innerHTML = `<strong>${escapeHtml(item.name)}</strong><span>${escapeHtml(item.source)}${item.selected ? " · selected" : ""}</span>`;
    this.tooltip.classList.remove("hidden");
    this.positionTooltip(event);
  }

  positionTooltip(event) {
    if (this.tooltip.classList.contains("hidden")) return;
    const rect = this.container.getBoundingClientRect();
    const rawX = (event.clientX || rect.left + rect.width / 2) - rect.left + 13;
    const rawY = (event.clientY || rect.top + rect.height / 2) - rect.top + 13;
    const width = this.tooltip.offsetWidth || 220;
    const height = this.tooltip.offsetHeight || 54;
    this.tooltip.style.left = `${clamp(rawX, 8, rect.width - width - 8)}px`;
    this.tooltip.style.top = `${clamp(rawY, 8, rect.height - height - 8)}px`;
  }

  hideTooltip() {
    this.tooltip.classList.add("hidden");
  }

  updateCoordinate(geo, suffix = "") {
    if (!geo) {
      this.coordinate.textContent = suffix || "Move over the map to inspect coordinates";
      return;
    }
    const base = `${formatHemisphere(geo.lat, "N", "S")} · ${formatHemisphere(geo.lon, "E", "W")}`;
    this.coordinate.textContent = suffix ? `${base} · ${suffix}` : base;
  }

  geoAtEvent(event) {
    const rect = this.container.getBoundingClientRect();
    const world = this.worldAtScreen(event.clientX - rect.left, event.clientY - rect.top);
    return unproject(world.x, world.y);
  }

  worldAtScreen(x, y) {
    const rect = this.container.getBoundingClientRect();
    const scale = this.worldScale();
    return {
      x: this.center.x + (x - rect.width / 2) / scale,
      y: clamp(this.center.y + (y - rect.height / 2) / scale, 0, 1),
    };
  }

  screenAtGeo(lon, lat) {
    const world = project(lon, lat);
    const scale = this.worldScale();
    let dx = world.x - this.center.x;
    dx -= Math.round(dx);
    return {
      x: this.container.clientWidth / 2 + dx * scale,
      y: this.container.clientHeight / 2 + (world.y - this.center.y) * scale,
    };
  }

  worldScale() {
    return 256 * (2 ** this.zoom);
  }

  queueRender() {
    window.cancelAnimationFrame(this.resizeFrame);
    this.resizeFrame = window.requestAnimationFrame(() => this.render());
  }
}

function project(lon, lat) {
  const latitude = clamp(Number(lat), -85.05112878, 85.05112878);
  const radians = latitude * Math.PI / 180;
  return {
    x: (Number(lon) + 180) / 360,
    y: (1 - Math.log(Math.tan(radians) + (1 / Math.cos(radians))) / Math.PI) / 2,
  };
}

function unproject(x, y) {
  const lon = wrapLongitude(x * 360 - 180);
  const n = Math.PI - 2 * Math.PI * clamp(y, 0, 1);
  const lat = 180 / Math.PI * Math.atan(Math.sinh(n));
  return { lon, lat: clamp(lat, -85.05112878, 85.05112878) };
}

function normaliseMatrix(value) {
  if (!value) return [];
  let matrix = value;
  while (Array.isArray(matrix) && matrix.length === 1 && Array.isArray(matrix[0]) && Array.isArray(matrix[0][0])) matrix = matrix[0];
  if (!Array.isArray(matrix)) return [];
  if (matrix.length === 2 && Array.isArray(matrix[0]) && Array.isArray(matrix[1]) && matrix[0].length > 2 && matrix[0].length === matrix[1].length) {
    matrix = matrix[0].map((lon, index) => [lon, matrix[1][index]]);
  }
  return matrix
    .filter(row => Array.isArray(row) && row.length >= 2)
    .map(row => [Number(row[0]), Number(row[1])])
    .filter(row => Number.isFinite(row[0]) && Number.isFinite(row[1]));
}

function asCollection(value) {
  if (!value) return [];
  return Array.isArray(value) ? value : [value];
}

function localTileUrl(layer, zoom, x, y) {
  const extension = layer === "labels" ? "png" : "jpg";
  return `map_tiles/${layer}/${zoom}/${x}/${y}.${extension}`;
}

function openPolygon(points) {
  const polygon = normaliseMatrix(points);
  if (polygon.length < 2) return polygon;
  const first = polygon[0], last = polygon[polygon.length - 1];
  if (first[0] === last[0] && first[1] === last[1]) polygon.pop();
  return polygon;
}

function closePolygon(points) {
  const polygon = removeRepeatedPoints(openPolygon(points));
  if (!polygon.length) return [];
  polygon.push([...polygon[0]]);
  return polygon;
}

function removeRepeatedPoints(points) {
  const result = [];
  normaliseMatrix(points).forEach(point => {
    const previous = result[result.length - 1];
    if (!previous || Math.abs(previous[0] - point[0]) + Math.abs(previous[1] - point[1]) > 1e-9) result.push(point);
  });
  return result;
}

function cloneMatrix(points) {
  return normaliseMatrix(points).map(point => [...point]);
}

function unwrapWorldXs(xs) {
  if (!xs.length) return [];
  const ordinary = [...xs];
  if (Math.max(...ordinary) - Math.min(...ordinary) >= 0.95) return ordinary;
  const wrapped = xs.map(x => x < 0.5 ? x + 1 : x);
  const span = values => Math.max(...values) - Math.min(...values);
  return span(wrapped) < span(ordinary) ? wrapped : ordinary;
}

function gridInterval(zoom) {
  if (zoom < 1.4) return 60;
  if (zoom < 2.4) return 30;
  if (zoom < 3.6) return 10;
  if (zoom < 5) return 5;
  if (zoom < 6.2) return 2;
  if (zoom < 7.5) return 1;
  if (zoom < 9) return 0.5;
  if (zoom < 10.5) return 0.2;
  return 0.1;
}

function formatHemisphere(value, positive, negative) {
  const direction = value >= 0 ? positive : negative;
  return `${Math.abs(value).toFixed(5)}°${direction}`;
}

function roundCoordinate(value) {
  return Math.round(Number(value) * 1e7) / 1e7;
}

function wrapLongitude(lon) {
  return ((Number(lon) + 180) % 360 + 360) % 360 - 180;
}

function wrap01(value) {
  return ((value % 1) + 1) % 1;
}

function clamp(value, minimum, maximum) {
  return Math.min(maximum, Math.max(minimum, value));
}

function svgNode(name, attributes = {}) {
  const node = document.createElementNS("http://www.w3.org/2000/svg", name);
  Object.entries(attributes).forEach(([key, value]) => node.setAttribute(key, value));
  return node;
}

function clearNode(node) {
  while (node.firstChild) node.removeChild(node.firstChild);
}

function escapeHtml(value) {
  const node = document.createElement("div");
  node.textContent = String(value ?? "");
  return node.innerHTML;
}

window.PhaseMap = PhaseMap;
