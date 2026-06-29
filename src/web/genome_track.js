// ── Global state ──────────────────────────────────────────────────────────────
const NC_COLOR  = '#6e7681';
const RULER_H   = 28;
const FH        = 19;
const LEVEL_GAP = 24;
const ARROW_H   = 14;
const STOP_COL  = '#f85149';
const STOP_SET  = new Set(['TAA','TAG','TGA']);
const COMP      = {A:'T',T:'A',G:'C',C:'G',a:'t',t:'a',g:'c',c:'g'};

let lastState = null, localVS = 0, localVE = 0;
let localMode = false, localModeTimer = null;
let pyodide = null, pyReady = false, pyRendering = false;
let lastSvgData = null, debounceTimer = null;
let ws, reconnTimer, dragStart = null, dragVS0 = 0, dragVE0 = 0;

// ── Helpers ───────────────────────────────────────────────────────────────────
function fmt(n) {
  if (n >= 1e6) return (n/1e6).toFixed(2)+' Mb';
  if (n >= 1e3) return (n/1e3).toFixed(1)+' kb';
  return n+' bp';
}
function trunc(s,n)  { return s.length<=n ? s : s.slice(0,n-1)+'\u2026'; }
function esc(s)      { return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;'); }
function setFigStatus(m) { document.getElementById('fig-status').textContent = m; }
function niceTick(span, W) {
  const raw = span / Math.max(Math.floor(W/90),1);
  const mag = Math.pow(10, Math.floor(Math.log10(Math.max(raw,1))));
  for (const m of [1,2,5,10]) { if (mag*m >= raw) return mag*m; }
  return mag*10;
}
function dlBlob(blob,name) {
  const url=URL.createObjectURL(blob);
  Object.assign(document.createElement('a'),{href:url,download:name}).click();
  URL.revokeObjectURL(url);
}
function setFig(html) { document.getElementById('fig-output').innerHTML=html; }

// ── Level assignment ──────────────────────────────────────────────────────────
function assignLevels(feats, scale, vs) {
  const ends = [];
  return feats.map(f => {
    const x1 = (f.start-vs)*scale, x2 = (f.end-vs)*scale;
    let lv = 0;
    while (ends[lv] !== undefined && ends[lv] > x1-6) lv++;
    ends[lv] = x2+6;
    return {...f, x1, x2, lv};
  });
}

// ── Arrow polygon ─────────────────────────────────────────────────────────────
function arrowPts(x1,x2,cy,h,strand) {
  const w=x2-x1, ah=Math.min(h*.9,w*.25,10), t=cy-h/2, b=cy+h/2;
  if (strand==='+') {
    if (w<3) return `${x1},${t} ${x2},${cy} ${x1},${b}`;
    return `${x1},${t} ${x2-ah},${t} ${x2},${cy} ${x2-ah},${b} ${x1},${b}`;
  } else {
    if (w<3) return `${x2},${t} ${x1},${cy} ${x2},${b}`;
    return `${x2},${t} ${x1+ah},${t} ${x1},${cy} ${x1+ah},${b} ${x2},${b}`;
  }
}

// ── Six-frame helpers ─────────────────────────────────────────────────────────
function rcSeq(s) {
  let r=''; for (let i=s.length-1;i>=0;i--) r += COMP[s[i]]||'N'; return r;
}
function stopPositions(seq, frame) {
  const out=[], u=seq.toUpperCase();
  for (let i=frame; i+2<u.length; i+=3) { if (STOP_SET.has(u.slice(i,i+3))) out.push(i); }
  return out;
}

// ── SVG render ────────────────────────────────────────────────────────────────
function drawSVG(state, vs, ve) {
  const svgEl = document.getElementById('svg');
  const panel = document.getElementById('live-panel');
  const W = panel.clientWidth  || window.innerWidth;
  const H = panel.clientHeight || 300;
  const {genome_size, genome_name, features, sequence='', seq_start: seqOff=0} = state;
  const span  = Math.max(ve-vs,1), scale = W/span;
  const sixFrame = document.getElementById('opt-sixframe').checked;
  const covStyle = document.getElementById('cov-style-svg').value;
  const hasCovData = state.coverage_plus && state.coverage_plus.length > 0;
  const covH = (covStyle !== 'none' && hasCovData)
    ? parseInt(document.getElementById('cov-height').value) : 0;

  document.getElementById('genome-name').textContent = genome_name || 'genome';
  document.getElementById('position').textContent =
    `${fmt(vs)} \u2013 ${fmt(ve)}  (${fmt(span)})`;

  const out = [];

  // Ruler
  const tick=niceTick(span,W), first=Math.ceil(vs/tick)*tick;
  out.push(`<line x1="0" y1="${RULER_H}" x2="${W}" y2="${RULER_H}" stroke="#30363d" stroke-width="1"/>`);
  for (let pos=first; pos<=ve; pos+=tick) {
    const x=(pos-vs)*scale;
    out.push(`<line x1="${x|0}" y1="${RULER_H-7}" x2="${x|0}" y2="${RULER_H}" stroke="#484f58" stroke-width="1"/>`);
    out.push(`<text x="${(x+3)|0}" y="${RULER_H-10}" font-size="10" fill="#8b949e" font-family="monospace">${fmt(pos)}</text>`);
  }

  const yGeneTop = RULER_H + covH;
  if (covH > 0) drawCoverageStrand(state, vs, ve, W, RULER_H, covH, true, out);
  else if (covStyle !== 'none' && !hasCovData) {
    out.push(`<rect x="0" y="${RULER_H}" width="${W}" height="20" fill="#0a0e18"/>`);
    out.push(`<text x="${W/2|0}" y="${RULER_H+13}" text-anchor="middle" font-size="10" fill="#3d4a5e" font-family="monospace">coverage: no BAM loaded (--bam reads.bam)</text>`);
  }

  let geneBottom;
  if (sixFrame) {
    geneBottom = drawSixFrame(features, sequence, seqOff, vs, ve, W, scale, out, yGeneTop);
  } else {
    geneBottom = drawSimple(features, vs, ve, W, scale, out, yGeneTop);
  }

  if (covH > 0) drawCoverageStrand(state, vs, ve, W, geneBottom, covH, false, out);

  // Minimap
  if (genome_size > 0) {
    const mw=Math.min(W-20,200), mh=4, mx=W-mw-8, my=H-12;
    out.push(`<rect x="${mx}" y="${my}" width="${mw}" height="${mh}" fill="#21262d" rx="2"/>`);
    const vx1=mx+(vs/genome_size)*mw, vx2=mx+(ve/genome_size)*mw;
    out.push(`<rect x="${vx1|0}" y="${my}" width="${Math.max(2,(vx2-vx1))|0}" height="${mh}" fill="#58a6ff" rx="1"/>`);
  }

  svgEl.innerHTML = out.join('');
}

// ── Simple +/- track ──────────────────────────────────────────────────────────
function drawSimple(features, vs, ve, W, scale, out, yTop) {
  const SGAP = 8;
  const plus  = assignLevels(features.filter(f=>f.strand==='+'), scale, vs);
  const minus = assignLevels(features.filter(f=>f.strand==='-'), scale, vs);
  const nPlus  = plus.length  ? Math.max(0,...plus.map(f=>f.lv))+1 : 0;
  const nMinus = minus.length ? Math.max(0,...minus.map(f=>f.lv))+1 : 0;
  const divY   = yTop + nPlus*LEVEL_GAP + SGAP;
  out.push(`<line x1="0" y1="${divY}" x2="${W}" y2="${divY}" stroke="#21262d" stroke-width="1"/>`);
  const drawG = (f, cy, strand) => {
    const col = f.noncoding ? NC_COLOR : f.color;
    const x1c=Math.max(0,f.x1), x2c=Math.min(W,f.x2);
    if (x2c<=x1c) return;
    const pts=arrowPts(x1c,x2c,cy,ARROW_H,strand), lw=x2c-x1c;
    out.push(`<g class="gene" data-idx="${f._idx}"><title>${esc(f.name)}</title><polygon points="${pts}" fill="${col}" fill-opacity="0.88"/>`);
    if (lw>24) {
      const lbl=trunc(f.name,Math.floor(lw/6.5));
      const ly = strand==='+'? (cy-ARROW_H/2-3):(cy+ARROW_H/2+11);
      out.push(`<text x="${((x1c+x2c)/2)|0}" y="${ly|0}" text-anchor="middle" font-size="10" fill="#c9d1d9" font-family="monospace">${esc(lbl)}</text>`);
    }
    out.push('</g>');
  };
  plus.forEach(f  => drawG(f, divY-SGAP-ARROW_H/2-f.lv*LEVEL_GAP, '+'));
  minus.forEach(f => drawG(f, divY+SGAP+ARROW_H/2+f.lv*LEVEL_GAP, '-'));
  return nMinus > 0
    ? (divY + SGAP + ARROW_H + (nMinus-1)*LEVEL_GAP + 14)|0
    : divY + SGAP*2;
}

// ── Six-frame translation track ───────────────────────────────────────────────
function drawSixFrame(features, seq, seqOff, vs, ve, W, scale, out, yTop) {
  const rc      = rcSeq(seq);
  const seqLen  = seq.length;
  const showSC  = document.getElementById('opt-stopcodons').checked;
  const fwdStops = [[],[],[]];
  const revStops = [[],[],[]];
  if (showSC && seq) {
    for (let rcFr = 0; rcFr < 3; rcFr++) {
      for (const p of stopPositions(seq, rcFr)) {
        const gp = seqOff + p;
        fwdStops[gp % 3].push(gp);
      }
      for (const p of stopPositions(rc, rcFr)) {
        const gEnd = seqOff + seqLen - p - 1;
        revStops[gEnd % 3].push(seqOff + seqLen - p - 2);
      }
    }
  }

  const divY = yTop + FH*3;
  const ncYp = divY + FH*3 + 3;
  const ncYm = ncYp + FH + 2;

  for (let fr=0;fr<3;fr++) {
    const bg = fr%2 ? '#0f1520' : '#0d1117';
    out.push(`<rect x="0" y="${yTop+FH*fr}" width="${W}" height="${FH}" fill="${bg}"/>`);
    out.push(`<rect x="0" y="${divY+FH*fr}" width="${W}" height="${FH}" fill="${bg}"/>`);
  }

  out.push(`<line x1="0" y1="${divY}" x2="${W}" y2="${divY}" stroke="#30363d" stroke-width="1"/>`);

  for (let fr=0;fr<3;fr++) {
    const y1=yTop+FH*fr+2, y2=yTop+FH*(fr+1)-2;
    for (const gp of fwdStops[fr]) {
      const x=(gp-vs)*scale;
      if (x<-1||x>W+1) continue;
      out.push(`<line x1="${x|0}" y1="${y1}" x2="${x|0}" y2="${y2}" stroke="${STOP_COL}" stroke-width="1" opacity="0.6"/>`);
    }
    out.push(`<text x="3" y="${(yTop+FH*fr+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">+${fr+1}</text>`);
  }
  for (let fr=0;fr<3;fr++) {
    const y1=divY+FH*fr+2, y2=divY+FH*(fr+1)-2;
    for (const gp of revStops[fr]) {
      const x=(gp-vs)*scale;
      if (x<-1||x>W+1) continue;
      out.push(`<line x1="${x|0}" y1="${y1}" x2="${x|0}" y2="${y2}" stroke="${STOP_COL}" stroke-width="1" opacity="0.6"/>`);
    }
    out.push(`<text x="3" y="${(divY+FH*fr+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">-${fr+1}</text>`);
  }

  const gh = FH-6;
  for (const f of features) {
    if (f.noncoding) continue;
    const x1=(f.start-vs)*scale, x2=(f.end-vs)*scale;
    const x1c=Math.max(0,x1), x2c=Math.min(W,x2);
    if (x2c<=x1c) continue;
    const col = f.color;
    const fr  = f.strand==='+' ? (f.start-1)%3 : (f.end-1)%3;
    const cy  = f.strand==='+' ? yTop+FH*fr+FH/2 : divY+FH*fr+FH/2;
    const pts = arrowPts(x1c, x2c, cy, gh, f.strand);
    const lw  = x2c-x1c;
    out.push(`<g class="gene" data-idx="${f._idx}"><title>${esc(f.name)}</title><polygon points="${pts}" fill="${col}" fill-opacity="0.9"/>`);
    if (lw>18) {
      const lbl = trunc(f.name, Math.floor(lw/5.5));
      out.push(`<text x="${((x1c+x2c)/2)|0}" y="${(cy+3)|0}" text-anchor="middle" font-size="9" fill="#fff" fill-opacity="0.85" font-family="monospace">${esc(lbl)}</text>`);
    }
    out.push('</g>');
  }

  out.push(`<line x1="0" y1="${ncYp-2}" x2="${W}" y2="${ncYp-2}" stroke="#21262d" stroke-width="1"/>`);
  out.push(`<text x="3" y="${(ncYp+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">nc+</text>`);
  out.push(`<text x="3" y="${(ncYm+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">nc-</text>`);
  for (const f of features) {
    if (!f.noncoding) continue;
    const x1=(f.start-vs)*scale, x2=(f.end-vs)*scale;
    const x1c=Math.max(0,x1), x2c=Math.min(W,x2);
    if (x2c<=x1c) continue;
    const cy = f.strand==='+' ? ncYp+FH/2 : ncYm+FH/2;
    const pts = arrowPts(x1c, x2c, cy, FH-6, f.strand);
    out.push(`<g class="gene" data-idx="${f._idx}"><title>${esc(f.name)}</title><polygon points="${pts}" fill="${NC_COLOR}" fill-opacity="0.7"/></g>`);
  }
  return ncYm + FH + 4;
}
