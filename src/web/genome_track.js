// ── Global state ──────────────────────────────────────────────────────────────
const NC_COLOR  = '#6e7681';
const RULER_H   = 28;
const FH        = 19;
const LEVEL_GAP = 24;
const ARROW_H   = 14;
const STOP_COL  = '#f85149';
const COMP      = {A:'T',T:'A',G:'C',C:'G',a:'t',t:'a',g:'c',c:'g'};

// ── Genetic code tables ────────────────────────────────────────────────────────
const GENETIC_CODES = {
  '1':  { stops: new Set(['TAA','TAG','TGA']) },
  '2':  { stops: new Set(['TAA','TAG','AGA','AGG']) },
  '4':  { stops: new Set(['TAA','TAG']) },
  '11': { stops: new Set(['TAA','TAG','TGA']) },
};
const AA_TABLE = {
  TTT:'F',TTC:'F',TTA:'L',TTG:'L',CTT:'L',CTC:'L',CTA:'L',CTG:'L',
  ATT:'I',ATC:'I',ATA:'I',ATG:'M',GTT:'V',GTC:'V',GTA:'V',GTG:'V',
  TCT:'S',TCC:'S',TCA:'S',TCG:'S',CCT:'P',CCC:'P',CCA:'P',CCG:'P',
  ACT:'T',ACC:'T',ACA:'T',ACG:'T',GCT:'A',GCC:'A',GCA:'A',GCG:'A',
  TAT:'Y',TAC:'Y',TAA:'*',TAG:'*',CAT:'H',CAC:'H',CAA:'Q',CAG:'Q',
  AAT:'N',AAC:'N',AAA:'K',AAG:'K',GAT:'D',GAC:'D',GAA:'E',GAG:'E',
  TGT:'C',TGC:'C',TGA:'*',TGG:'W',CGT:'R',CGC:'R',CGA:'R',CGG:'R',
  AGT:'S',AGC:'S',AGA:'R',AGG:'R',GGT:'G',GGC:'G',GGA:'G',GGG:'G',
};
function getStopSet() {
  const sel = document.getElementById('opt-gencode');
  return (GENETIC_CODES[(sel && sel.value) || '1'] || GENETIC_CODES['1']).stops;
}
function translateCodon(c) { return AA_TABLE[c.toUpperCase()] || '?'; }

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
function stopPositions(seq, frame, stopSet) {
  const out=[], u=seq.toUpperCase(), s=stopSet||getStopSet();
  for (let i=frame; i+2<u.length; i+=3) { if (s.has(u.slice(i,i+3))) out.push(i); }
  return out;
}

// ── Arrow style (set by fancy_track.js; fallback for load order) ─────────────
if (!window.getFancyMode) window.getFancyMode = function() { return false; };

// ── SVG render ────────────────────────────────────────────────────────────────
function drawSVG(state, vs, ve) {
  const svgEl = document.getElementById('svg');
  const panel = document.getElementById('live-panel');
  const W = panel.clientWidth  || window.innerWidth;
  const H = panel.clientHeight || 300;
  const {genome_size, genome_name, sequence='', seq_start: seqOff=0} = state;
  // Merge blast hits (yellow) into the feature list — overlaid in same gene track row
  const blastFeats = (state.blast_features || []).map(bf => ({...bf, color: '#e6b800'}));
  const features = [...(state.features || []), ...blastFeats];
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
  window._fancyYGeneTop = yGeneTop;  // read by fancy_track.js for tile positioning
  if (covH > 0) drawCoverageStrand(state, vs, ve, W, RULER_H, covH, true, out);
  else if (covStyle !== 'none' && !hasCovData) {
    out.push(`<rect x="0" y="${RULER_H}" width="${W}" height="20" fill="#0a0e18"/>`);
    out.push(`<text x="${W/2|0}" y="${RULER_H+13}" text-anchor="middle" font-size="10" fill="#3d4a5e" font-family="monospace">coverage: no BAM loaded (--bam reads.bam)</text>`);
  }

  const hitOnly = window.getFancyMode();
  let geneBottom;
  if (sixFrame) {
    geneBottom = drawSixFrame(features, sequence, seqOff, vs, ve, W, scale, out, yGeneTop, hitOnly);
  } else {
    geneBottom = drawSimple(features, vs, ve, W, scale, out, yGeneTop, hitOnly);
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
function drawSimple(features, vs, ve, W, scale, out, yTop, hitOnly=false) {
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
    if (hitOnly) {
      // Transparent hit target — tiles provide the visible arrow
      out.push(`<g class="gene" data-idx="${f._idx}"><title>${esc(f.name)}</title><polygon points="${pts}" fill="transparent" stroke="none" pointer-events="all"/></g>`);
    } else {
      out.push(`<g class="gene" data-idx="${f._idx}"><title>${esc(f.name)}</title><polygon points="${pts}" fill="${col}" fill-opacity="0.88"/>`);
      if (lw>24) {
        const lbl=trunc(f.name,Math.floor(lw/6.5));
        const ly = strand==='+'? (cy-ARROW_H/2-3):(cy+ARROW_H/2+11);
        out.push(`<text x="${((x1c+x2c)/2)|0}" y="${ly|0}" text-anchor="middle" font-size="10" fill="#c9d1d9" font-family="monospace">${esc(lbl)}</text>`);
      }
      out.push('</g>');
    }
  };
  plus.forEach(f  => drawG(f, divY-SGAP-ARROW_H/2-f.lv*LEVEL_GAP, '+'));
  minus.forEach(f => drawG(f, divY+SGAP+ARROW_H/2+f.lv*LEVEL_GAP, '-'));
  return nMinus > 0
    ? (divY + SGAP + ARROW_H + (nMinus-1)*LEVEL_GAP + 14)|0
    : divY + SGAP*2;
}

// ── Nucleotide row ────────────────────────────────────────────────────────────
function drawNucleotideRow(seq, seqOff, vs, W, scale, out, yTop, rowH, isComplement) {
  const BASE_COL = {A:'#2ea043',T:'#da3633',G:'#388bfd',C:'#e3b341'};
  out.push(`<rect x="0" y="${yTop}" width="${W}" height="${rowH}" fill="#0a0d14"/>`);
  const label = isComplement ? '-nt' : '+nt';
  out.push(`<text x="3" y="${(yTop+rowH/2+3)|0}" font-size="7" fill="#484f58" font-family="monospace">${label}</text>`);
  if (!seq || scale < 1.5) return;
  const u = seq.toUpperCase();
  for (let i = 0; i < u.length; i++) {
    const gpos = seqOff + i;
    const x = (gpos - vs) * scale;
    if (x + scale < 0 || x > W) continue;
    const rawBase = u[i];
    const base = isComplement ? (COMP[rawBase] || 'N') : rawBase;
    const col = BASE_COL[base] || '#484f58';
    const bx  = Math.max(0, x|0);
    const bx2 = Math.min(W, (x + scale)|0);
    if (bx2 <= bx) continue;
    out.push(`<rect x="${bx}" y="${yTop+1}" width="${bx2-bx}" height="${rowH-2}" fill="${col}" opacity="0.4"/>`);
    if (scale >= 9) {
      const fs = Math.min(rowH - 2, (scale * 0.65)|0);
      if (fs >= 5)
        out.push(`<text x="${(x+scale/2)|0}" y="${(yTop+(rowH+fs)*0.5)|0}" text-anchor="middle" font-size="${fs}" fill="${col}" font-family="monospace">${base}</text>`);
    }
  }
}

// ── Six-frame translation track ───────────────────────────────────────────────
// Layout (top → bottom):
//   +fr0 AA/stop | +fr0 genes | +fr1 AA/stop | +fr1 genes | +fr2 AA/stop | +fr2 genes
//   [+nt] | separator | [-nt]
//   -fr0 AA/stop | -fr0 genes | -fr1 AA/stop | -fr1 genes | -fr2 AA/stop | -fr2 genes
//   nc+ | nc-
function drawSixFrame(features, seq, seqOff, vs, ve, W, scale, out, yTop, hitOnly=false) {
  const GH     = 13;   // gene-arrow row height per frame
  const NH     = 13;   // nucleotide row height
  const SEP    = 4;    // strand separator height
  const showNT = !!(document.getElementById('opt-nucleotides') || {checked:false}).checked;

  // ── Row Y positions ────────────────────────────────────────────────────────
  let y = yTop;
  const fwdAaY = [], fwdGnY = [];
  for (let fr = 0; fr < 3; fr++) {
    fwdAaY[fr] = y; y += FH;
    fwdGnY[fr] = y; y += GH;
  }
  let ntPlusY = -1;
  if (showNT) { ntPlusY = y; y += NH; }
  const sepY = y; y += SEP;
  let ntMinusY = -1;
  if (showNT) { ntMinusY = y; y += NH; }
  const revAaY = [], revGnY = [];
  for (let fr = 0; fr < 3; fr++) {
    revAaY[fr] = y; y += FH;
    revGnY[fr] = y; y += GH;
  }
  const ncYp = y; y += FH;
  const ncYm = y; y += FH;
  const totalBottom = y;

  // ── hitOnly: transparent click targets covering full height ────────────────
  if (hitOnly) {
    const tileH = totalBottom - yTop;
    for (const f of features) {
      const x1c = Math.max(0, (f.start - vs) * scale);
      const x2c = Math.min(W,  (f.end   - vs) * scale);
      if (x2c <= x1c) continue;
      out.push(`<g class="gene" data-idx="${f._idx}"><title>${esc(f.name)}</title>`
        + `<rect x="${x1c|0}" y="${yTop}" width="${(x2c-x1c)|0}" height="${tileH}" `
        + `fill="transparent" stroke="none" pointer-events="all"/></g>`);
    }
    return totalBottom;
  }

  // ── Collect stop codons ────────────────────────────────────────────────────
  const rc      = rcSeq(seq);
  const seqLen  = seq.length;
  const showSC  = document.getElementById('opt-stopcodons').checked;
  const stopSet = getStopSet();
  // fwdStops[fr] = [gp_left, …]  where gp_left is the first (leftmost) base of the codon
  // revStops[fr] = [gp_left, …]  same convention; frame keyed by rightmost base % 3
  const fwdStops = [[],[],[]];
  const revStops = [[],[],[]];
  if (showSC) {
    const prePlus  = lastState && lastState.stop_codons_plus  && lastState.stop_codons_plus.length  === 3;
    const preMinus = lastState && lastState.stop_codons_minus && lastState.stop_codons_minus.length === 3;
    if (prePlus && preMinus) {
      // Server sends 1-based positions, pre-bucketed by frame.
      // fwdStops[fr] = 0-based gpLeft; revStops[fr] = 0-based gLeft (leftmost base of codon).
      for (let fr = 0; fr < 3; fr++) {
        for (const pos1 of lastState.stop_codons_plus[fr])  fwdStops[fr].push(pos1 - 1);
        for (const pos1 of lastState.stop_codons_minus[fr]) revStops[fr].push(pos1 - 1);
      }
    } else if (seq) {
      // Fallback: scan sequence window (limited to ±500kb)
      for (let rcFr = 0; rcFr < 3; rcFr++) {
        for (const p of stopPositions(seq, rcFr, stopSet)) {
          const gpLeft = seqOff + p;
          fwdStops[gpLeft % 3].push(gpLeft);
        }
        for (const p of stopPositions(rc, rcFr, stopSet)) {
          const gEnd  = seqOff + seqLen - 1 - p;
          const gLeft = gEnd - 2;
          revStops[gEnd % 3].push(gLeft);
        }
      }
    }
  }

  // ── Background rects ───────────────────────────────────────────────────────
  for (let fr = 0; fr < 3; fr++) {
    const bg  = fr % 2 ? '#0f1520' : '#0d1117';
    const bgG = fr % 2 ? '#0c1019' : '#0b0f16';
    out.push(`<rect x="0" y="${fwdAaY[fr]}" width="${W}" height="${FH}"  fill="${bg}"/>`);
    out.push(`<rect x="0" y="${fwdGnY[fr]}" width="${W}" height="${GH}"  fill="${bgG}"/>`);
    out.push(`<rect x="0" y="${revAaY[fr]}" width="${W}" height="${FH}"  fill="${bg}"/>`);
    out.push(`<rect x="0" y="${revGnY[fr]}" width="${W}" height="${GH}"  fill="${bgG}"/>`);
  }
  out.push(`<rect x="0" y="${ncYp}" width="${W}" height="${FH}" fill="#0d1117"/>`);
  out.push(`<rect x="0" y="${ncYm}" width="${W}" height="${FH}" fill="#0d1117"/>`);

  // Strand separator line
  out.push(`<line x1="0" y1="${sepY+1}" x2="${W}" y2="${sepY+1}" stroke="#30363d" stroke-width="1"/>`);
  // nc separator
  out.push(`<line x1="0" y1="${ncYp-1}" x2="${W}" y2="${ncYp-1}" stroke="#21262d" stroke-width="1"/>`);

  // ── AA translation text ────────────────────────────────────────────────────
  const showAA = seq && scale >= 3;
  const aaFont = Math.min(11, Math.max(7, (scale * 2.2)|0));
  if (showAA) {
    const u    = seq.toUpperCase();
    const u_rc = rc.toUpperCase();
    // Forward: codons whose leftmost base ≡ fr (mod 3)
    for (let fr = 0; fr < 3; fr++) {
      const cy     = (fwdAaY[fr] + FH/2 + 3)|0;
      const pStart = ((fr - seqOff % 3 + 3) % 3);
      for (let p = pStart; p + 2 < seqLen; p += 3) {
        const gp = seqOff + p;
        const cx = (gp + 1.5 - vs) * scale;
        if (cx < -20 || cx > W + 20) continue;
        const aa = translateCodon(u.slice(p, p + 3));
        if (!aa || aa === '*') continue;
        out.push(`<text x="${cx|0}" y="${cy}" text-anchor="middle" font-size="${aaFont}" fill="#4a5568" font-family="monospace">${aa}</text>`);
      }
    }
    // Reverse: codons keyed by rightmost-base frame
    for (let rcFr = 0; rcFr < 3; rcFr++) {
      for (let p = rcFr; p + 2 < seqLen; p += 3) {
        const gEnd  = seqOff + seqLen - 1 - p;
        const gLeft = gEnd - 2;
        const cx    = (gLeft + 1.5 - vs) * scale;
        if (cx < -20 || cx > W + 20) continue;
        const fr = gEnd % 3;
        const cy = (revAaY[fr] + FH/2 + 3)|0;
        const aa = translateCodon(u_rc.slice(p, p + 3));
        if (!aa || aa === '*') continue;
        out.push(`<text x="${cx|0}" y="${cy}" text-anchor="middle" font-size="${aaFont}" fill="#4a5568" font-family="monospace">${aa}</text>`);
      }
    }
  }

  // ── Stop codons: centred box + asterisk ────────────────────────────────────
  // Box is exactly 3 bp wide, centred on codon midpoint (gLeft + 1.5 bp).
  const stopW = Math.max(2, scale * 3);
  for (let fr = 0; fr < 3; fr++) {
    const y1f = fwdAaY[fr] + 2, y2f = fwdAaY[fr] + FH - 2;
    for (const gpLeft of fwdStops[fr]) {
      const cx  = (gpLeft + 1.5 - vs) * scale;
      if (cx + stopW/2 < 0 || cx - stopW/2 > W) continue;
      const rx  = Math.max(0, (cx - stopW/2)|0);
      const rw  = Math.min(stopW, W - rx);
      out.push(`<rect x="${rx}" y="${y1f}" width="${rw}" height="${y2f-y1f}" fill="${STOP_COL}" opacity="0.28" rx="1"/>`);
      if (showAA)
        out.push(`<text x="${cx|0}" y="${(fwdAaY[fr]+FH/2+3)|0}" text-anchor="middle" font-size="${aaFont}" fill="${STOP_COL}" font-family="monospace">*</text>`);
    }
    out.push(`<text x="3" y="${(fwdAaY[fr]+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">+${fr+1}</text>`);
  }
  for (let fr = 0; fr < 3; fr++) {
    const y1r = revAaY[fr] + 2, y2r = revAaY[fr] + FH - 2;
    for (const gpLeft of revStops[fr]) {
      const cx  = (gpLeft + 1.5 - vs) * scale;
      if (cx + stopW/2 < 0 || cx - stopW/2 > W) continue;
      const rx  = Math.max(0, (cx - stopW/2)|0);
      const rw  = Math.min(stopW, W - rx);
      out.push(`<rect x="${rx}" y="${y1r}" width="${rw}" height="${y2r-y1r}" fill="${STOP_COL}" opacity="0.28" rx="1"/>`);
      if (showAA)
        out.push(`<text x="${cx|0}" y="${(revAaY[fr]+FH/2+3)|0}" text-anchor="middle" font-size="${aaFont}" fill="${STOP_COL}" font-family="monospace">*</text>`);
    }
    out.push(`<text x="3" y="${(revAaY[fr]+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">-${fr+1}</text>`);
  }

  // ── Gene arrows in their per-frame rows ────────────────────────────────────
  const gh = GH - 4;
  for (const f of features) {
    if (f.noncoding) continue;
    const x1c = Math.max(0, (f.start - vs) * scale);
    const x2c = Math.min(W,  (f.end   - vs) * scale);
    if (x2c <= x1c) continue;
    const fr  = f.strand === '+' ? (f.start - 1) % 3 : (f.end - 1) % 3;
    const gnY = f.strand === '+' ? fwdGnY[fr] : revGnY[fr];
    const cy  = gnY + GH / 2;
    const pts = arrowPts(x1c, x2c, cy, gh, f.strand);
    const lw  = x2c - x1c;
    out.push(`<g class="gene" data-idx="${f._idx}"><title>${esc(f.name)}</title><polygon points="${pts}" fill="${f.color}" fill-opacity="0.9"/>`);
    if (lw > 18) {
      out.push(`<text x="${((x1c+x2c)/2)|0}" y="${(cy+3)|0}" text-anchor="middle" font-size="8" fill="#fff" fill-opacity="0.85" font-family="monospace">${esc(trunc(f.name, Math.floor(lw/5.5)))}</text>`);
    }
    out.push('</g>');
  }

  // ── Nucleotide tracks ──────────────────────────────────────────────────────
  if (showNT) {
    drawNucleotideRow(seq, seqOff, vs, W, scale, out, ntPlusY,  NH, false);
    drawNucleotideRow(seq, seqOff, vs, W, scale, out, ntMinusY, NH, true);
  }

  // ── Non-coding rows ────────────────────────────────────────────────────────
  out.push(`<text x="3" y="${(ncYp+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">nc+</text>`);
  out.push(`<text x="3" y="${(ncYm+FH/2+3)|0}" font-size="8" fill="#484f58" font-family="monospace">nc-</text>`);
  for (const f of features) {
    if (!f.noncoding) continue;
    const x1c = Math.max(0, (f.start - vs) * scale);
    const x2c = Math.min(W,  (f.end   - vs) * scale);
    if (x2c <= x1c) continue;
    const cy  = (f.strand === '+' ? ncYp : ncYm) + FH / 2;
    const pts = arrowPts(x1c, x2c, cy, FH - 6, f.strand);
    out.push(`<g class="gene" data-idx="${f._idx}"><title>${esc(f.name)}</title><polygon points="${pts}" fill="${NC_COLOR}" fill-opacity="0.7"/></g>`);
  }

  return totalBottom;
}
