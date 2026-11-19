"use strict";

/* Geographic screening cells are drawn row by row in Mercator, so a WGS84
 * raster is never stretched as though latitude were a Mercator coordinate. */
class PhaseQualityMap extends PhaseMap {
  constructor(container, options = {}) {
    super(container, options);
    this.qualityCanvas = document.createElement("canvas");
    Object.assign(this.qualityCanvas.style, { position: "absolute", inset: "0", pointerEvents: "none" });
    this.container.insertBefore(this.qualityCanvas, this.labelLayer);
    this.qualityData = null;
    this.qualityRaster = null;
    this.qualityOpacity = 0.55;
    this.qualityLayer = "suitability";
  }

  fitBounds(points, animate = true) {
    const projected = normaliseMatrix(points).map(p => project(p[0], p[1]));
    if (!projected.length) return;
    const xs = unwrapWorldXs(projected.map(p => p.x)), ys = projected.map(p => p.y);
    const minX = Math.min(...xs), maxX = Math.max(...xs), minY = Math.min(...ys), maxY = Math.max(...ys);
    const width = Math.max(320, this.container.clientWidth || 800), height = Math.max(260, this.container.clientHeight || 480);
    const padding = Math.min(82, Math.min(width, height)*0.12);
    const scale = Math.min((width-2*padding)/Math.max(maxX-minX, 1e-7), (height-2*padding)/Math.max(maxY-minY, 1e-7));
    this.center = {x:wrap01((minX+maxX)/2), y:clamp((minY+maxY)/2, 0.002, 0.998)};
    this.zoom = clamp(Math.log2(scale/256), 0.45, this.maxZoom); this.render();
  }

  setQuality(data, layer = "suitability", opacity = 0.55) {
    this.qualityData = data?.available ? data : null;
    this.qualityLayer = layer;
    this.qualityOpacity = opacity;
    this.qualityRaster = null;
    if (this.qualityData && layer !== "off") {
      const g = data.grid, values = data.layers[layer]?.values;
      if (values) {
        const raster = document.createElement("canvas");
        raster.width = g.columns; raster.height = g.rows;
        const ctx = raster.getContext("2d"), pixels = ctx.createImageData(g.columns, g.rows);
        const flat = values.flat(Infinity);
        flat.forEach((v, i) => {
          if (v === null || !Number.isFinite(v)) return;
          const rgb = qualityColour(layer, v);
          pixels.data.set([...rgb, 255], 4*i);
        });
        ctx.putImageData(pixels, 0, 0); this.qualityRaster = raster;
      }
    }
    this.render();
  }

  render() {
    super.render();
    if (!this.qualityCanvas) return;
    const width = this.container.clientWidth, height = this.container.clientHeight;
    this.qualityCanvas.width = width; this.qualityCanvas.height = height;
    const ctx = this.qualityCanvas.getContext("2d");
    if (!this.qualityRaster || !this.qualityData) return;
    ctx.globalAlpha = this.qualityOpacity; ctx.imageSmoothingEnabled = false;
    const g = this.qualityData.grid, [w, e, s, n] = g.bounds;
    for (let row = 0; row < g.rows; row++) {
      const a = this.screenAtGeo(w, n-(n-s)*row/g.rows);
      const b = this.screenAtGeo(e, n-(n-s)*(row+1)/g.rows);
      ctx.drawImage(this.qualityRaster, 0, row, g.columns, 1, a.x, a.y, b.x-a.x, b.y-a.y+0.25);
    }
  }

  updateCoordinate(geo, suffix = "") {
    super.updateCoordinate(geo, suffix);
    if (!geo || suffix || !this.qualityData || this.qualityLayer === "off") return;
    const d = this.qualityData, g = d.grid, [w, e, s, n] = g.bounds;
    if (geo.lon < w || geo.lon > e || geo.lat < s || geo.lat > n) return;
    const col = Math.min(g.columns-1, Math.floor((geo.lon-w)/(e-w)*g.columns));
    const row = Math.min(g.rows-1, Math.floor((n-geo.lat)/(n-s)*g.rows));
    const index = row*g.columns+col, value = d.layers[this.qualityLayer].values.flat(Infinity)[index];
    let label = "Not evaluated";
    if (value !== null && Number.isFinite(value)) {
      label = this.qualityLayer === "suitability" ? ["Not evaluated", "Favourable", "Intermediate", "Unfavourable"][value] : (this.qualityLayer === "coherence" ? value.toFixed(2) : `${(100*value).toFixed(0)}%`);
    }
    const reasons = d.reasonCodes.flat(Infinity)[index] || 0;
    const descriptions = qualityReasons(d.reasonBits, reasons);
    this.coordinate.textContent += ` · ${label}${descriptions.length ? ` · ${descriptions.join("; ")}` : ""}`;
  }
}

class PhaseQualityControls {
  constructor(container, map) {
    this.container = container; this.map = map; this.data = null;
    this.select = container.querySelector("[data-quality-layer]");
    this.opacity = container.querySelector("[data-quality-opacity]");
    this.legend = container.querySelector("[data-quality-legend]");
    this.message = container.querySelector("[data-quality-message]");
    this.select.addEventListener("change", () => this.apply());
    this.opacity.addEventListener("input", () => this.apply());
  }
  update(data) {
    this.data = data;
    this.select.disabled = !data?.available; this.opacity.disabled = !data?.available;
    this.message.textContent = data?.available ? `Master ${data.master} · ${data.acquisitionCount} acquisitions · ${data.interferogramCount} interferogram${data.interferogramCount === 1 ? "" : "s"}. ${data.message} Geometry: ${data.geometry.message}` : (data?.message || "Complete preprocessing or build a quality map for the existing export.");
    this.apply();
  }
  apply() {
    this.map.setQuality(this.data, this.select.value, Number(this.opacity.value)/100);
    const layer = this.select.value;
    this.legend.replaceChildren();
    if (!this.data?.available || layer === "off") return;
    const labels = layer === "suitability" ? ["Not evaluated", "Favourable", "Intermediate", "Unfavourable"] : (layer === "stability" ? ["0%", "≥5%"] : layer === "coherence" ? ["0", "1"] : ["0%", "100%"]);
    labels.forEach((text, i) => {
      const span = document.createElement("span"), dot = document.createElement("i");
      dot.style.background = `rgb(${qualityColour(layer, layer === "suitability" ? i : i)})`;
      span.append(dot, document.createTextNode(text)); this.legend.appendChild(span);
    });
    const description = this.data.layers[layer]?.description;
    this.legend.title = description || "";
  }
}

function qualityColour(layer, value) {
  const grey = [155, 164, 179], green = [31, 160, 110], amber = [236, 166, 36], red = [212, 65, 75];
  if (layer === "suitability") return [grey, green, amber, red][value] || grey;
  let t = Math.max(0, Math.min(1, value));
  if (layer === "geometry") t = 1-t;
  // Stable-pixel fractions are often small: keep the declared favourable
  // screening threshold visible without changing the underlying exported data.
  if (layer === "stability") t = Math.min(1, t/0.05);
  const a = t < 0.5 ? red : amber, b = t < 0.5 ? amber : green, f = t < 0.5 ? 2*t : 2*t-1;
  return a.map((v, i) => Math.round(v+(b[i]-v)*f));
}

function qualityReasons(bits, code) {
  const items = Array.isArray(bits) ? bits : (bits ? [bits] : []);
  return items.filter(item => code & Number(item.bit)).map(item => item.message);
}
