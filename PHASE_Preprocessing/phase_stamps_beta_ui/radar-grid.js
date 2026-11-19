"use strict";

const radarStageNames = ["", "Rejected during selection", "Rejected during weeding", "Retained after weeding"];
const radarStageColours = ["", "#a7acb4", "#f29936", "#df408f"];

class PhaseRadarMap extends PhaseMap {
  constructor(container, options) {
    super(container, options);
    this.viewportChanged = options.viewportChanged;
    this.canvas = document.createElement("canvas");
    Object.assign(this.canvas.style, {position:"absolute", inset:"0", zIndex:"2", pointerEvents:"none"});
    container.appendChild(this.canvas);
    this.frame = null; this.colour = "stage"; this.showGrid = true; this.hitPoints = [];
    container.addEventListener("pointermove", event => this.inspectPoint(event));
  }
  fitBounds(points) { PhaseQualityMap.prototype.fitBounds.call(this, points); }
  bounds() {
    const a=this.worldAtScreen(0,0), b=this.worldAtScreen(this.container.clientWidth,this.container.clientHeight);
    const nw=unproject(a.x,a.y), se=unproject(b.x,b.y);
    return [Math.max(-180,nw.lon),Math.min(180,se.lon),Math.max(-85,se.lat),Math.min(85,nw.lat)];
  }
  render() {
    super.render();
    if (!this.canvas) return;
    const w=this.container.clientWidth, h=this.container.clientHeight, d=Math.min(window.devicePixelRatio||1,2);
    this.canvas.width=w*d; this.canvas.height=h*d; this.canvas.style.width=`${w}px`;this.canvas.style.height=`${h}px`;
    const ctx=this.canvas.getContext("2d"); ctx.scale(d,d); this.hitPoints=[];
    if (this.frame) {
      if(this.showGrid) {
        ctx.strokeStyle="rgba(245,249,255,.65)";ctx.lineWidth=.65;
        for(const line of this.frame.lines||[]) {
          ctx.beginPath();let start=true;
          for(const pair of line) {
            if(!pair || pair[0]===null || pair[1]===null) {start=true;continue;}
            const p=this.screenAtGeo(pair[0],pair[1]);
            if(start) ctx.moveTo(p.x,p.y);else ctx.lineTo(p.x,p.y);start=false;
          }
          ctx.stroke();
        }
      }
      let rows=this.frame.points||[];if(rows.length && !Array.isArray(rows[0])) rows=[rows];
      for(const row of rows) {
        const p=this.screenAtGeo(row[0],row[1]);if(p.x<0||p.x>w||p.y<0||p.y>h) continue;
        if(this.colour==="da") {
          if(row[2]===null) ctx.fillStyle="#a7acb4";
          else {const t=Math.max(0,Math.min(1,row[2]/Math.max(this.frame.threshold,.001)));ctx.fillStyle=`hsl(${240-240*t},80%,55%)`;}
        } else ctx.fillStyle=radarStageColours[row[3]]||"#a7acb4";
        ctx.beginPath();ctx.arc(p.x,p.y,2.8,0,Math.PI*2);ctx.fill();
        this.hitPoints.push({x:p.x,y:p.y,row});
      }
      this.viewportChanged?.();
    }
  }
  inspectPoint(event) {
    if(this.pan) return;
    const r=this.container.getBoundingClientRect(), x=event.clientX-r.left,y=event.clientY-r.top;
    let best=null, distance=64;
    for(const p of this.hitPoints) {const d=(p.x-x)**2+(p.y-y)**2;if(d<distance){best=p.row;distance=d;}}
    if(best) this.updateCoordinate({lon:best[0],lat:best[1]},
      `${radarStageNames[best[3]]} · D_A ${best[2]===null?"unavailable":Number(best[2]).toFixed(3)} · az ${best[4]}, range ${best[5]}`);
  }
}

class PhaseRadarGridView {
  constructor(sendEvent) {
    this.send=sendEvent;this.busy=false;this.ready=false;this.running=false;this.lastKey="";this.workDir="";this.exportBusy=false;
    this.map=new PhaseRadarMap(document.getElementById("radar-map"), {
      onTilesRequested:requests=>this.send("MapTilesRequested",{requests}),
      viewportChanged:()=>this.queue(),
    });
    this.el=id=>document.getElementById(id);
    this.el("grid-refresh").onclick=()=>this.load(true);
    this.el("grid-fit").onclick=()=>this.fit();
    this.el("grid-open").onclick=()=>this.send("OpenGridExports",{});
    this.el("grid-layer").onchange=()=>{this.map.colour=this.el("grid-layer").value;this.map.render();this.legend();};
    this.el("grid-stage").onchange=()=>{this.lastKey="";this.queue();};
    this.el("grid-lines").onchange=()=>{this.map.showGrid=this.el("grid-lines").checked;this.map.render();};
    this.el("grid-export").onclick=()=>{
      if(this.busy||this.running||!this.ready||this.exportBusy)return;
      this.exportBusy=true;this.controls();this.message("Exporting PNG/FIG maps and the complete candidate table…");
      this.send("ExportGridMap",{bounds:this.map.bounds()});
    };
    this.legend();
  }
  setState(state) {
    this.running=Boolean(state.running);
    if(this.workDir && (this.workDir!==state.workDir || this.running)) {
      this.ready=false;this.lastKey="";this.map.frame=null;this.map.render();
    }
    this.workDir=state.workDir;this.controls();
  }
  controls() {
    this.el("grid-refresh").disabled=this.running||this.busy||this.exportBusy;
    this.el("grid-export").disabled=this.running||this.busy||!this.ready||this.exportBusy;
  }
  message(text) {this.el("grid-message").textContent=text;}
  show() {this.map.render();if(!this.ready&&!this.busy&&!this.running)this.load(false);}
  load(refresh) {
    if(this.running||this.busy||this.exportBusy)return;
    this.ready=false;this.lastKey="";this.busy=true;this.controls();
    this.message("Loading radar geolocation and candidate selection history…");
    this.send("RefreshGridMap",{refresh,stage:Number(this.el("grid-stage").value)});
  }
  fit() {if(this.map.frame){const b=this.map.frame.extent;this.map.fitBounds([[b[0],b[2]],[b[1],b[3]]]);}}
  queue() {
    clearTimeout(this.timer);this.timer=setTimeout(()=>this.request(),260);
  }
  request() {
    if(!this.ready||this.busy||this.running||this.exportBusy||!document.getElementById("grid-page").classList.contains("active"))return;
    const bounds=this.map.bounds(),stage=Number(this.el("grid-stage").value);
    const key=JSON.stringify([bounds.map(x=>Math.round(x*1e7)/1e7),stage]);if(key===this.lastKey)return;
    this.lastKey=key;this.busy=true;this.controls();
    this.send("RefreshGridMap",{bounds,stage});
  }
  receive(data) {
    const initial=!this.ready;this.busy=false;
    if(!data?.available){this.ready=false;this.message(data?.message||"Radar-grid data unavailable.");this.controls();return;}
    this.ready=true;this.map.frame=data;
    if(initial)this.fit();else this.map.render();
    const sample=data.visibleCount>data.displayedCount?` Displaying ${data.displayedCount.toLocaleString()} sampled points; zoom in to inspect more.`:"";
    this.message(`${data.visibleCount.toLocaleString()} candidates in view. Grid lines every ${data.azimuthStep} azimuth / ${data.rangeStep} range cell(s). Geolocation: ${data.geolocationMode.replaceAll("_"," ")}.${sample}`);
    this.controls();this.queue();
  }
  exported(data) {this.exportBusy=false;this.message(data.message);this.controls();this.queue();}
  legend() {
    const node=this.el("grid-legend");node.replaceChildren();
    if(this.map.colour==="da") {node.textContent="D_A: blue = 0 → red = selection threshold or greater; grey = unavailable.";return;}
    [1,2,3].forEach(stage=>{const item=document.createElement("span"),dot=document.createElement("i");dot.style.background=radarStageColours[stage];item.append(dot,document.createTextNode(radarStageNames[stage]));node.append(item);});
  }
}
