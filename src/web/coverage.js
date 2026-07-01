// ── Coverage track ────────────────────────────────────────────────────────────
function gaussSmooth1D(arr, sigma) {
  const r=Math.ceil(3*sigma); const kern=[]; let sum=0;
  for (let i=-r;i<=r;i++) { const w=Math.exp(-0.5*(i/sigma)**2); kern.push(w); sum+=w; }
  kern.forEach((_,i,a)=>a[i]/=sum);
  return arr.map((_,ci)=>{
    let s=0;
    for (let j=0;j<kern.length;j++) { const idx=ci+j-r; if (idx>=0&&idx<arr.length) s+=arr[idx]*kern[j]; }
    return s;
  });
}

// Draw one strand's coverage track. isPlus=true → bars grow up; isPlus=false → bars grow down.
function drawCoverageStrand(state, vs, ve, W, yBase, trackH, isPlus, out) {
  const style = document.getElementById('cov-style-svg').value;
  if (style === 'none') return;
  const binSz  = state.coverage_bin_size  || 1000;
  const binOff = state.coverage_bin_start || 0;
  const data   = isPlus ? (state.coverage_plus||[]) : (state.coverage_minus||[]);
  if (!data.length) return;
  const span = Math.max(ve-vs,1), scaleX = W/span;
  const col  = isPlus ? '#58a6ff' : '#f78166';

  out.push(`<rect x="0" y="${yBase}" width="${W}" height="${trackH}" fill="#080d1a" opacity="0.9"/>`);

  const firstBin = Math.max(0, Math.floor((vs-binOff)/binSz));
  const lastBin  = Math.min(data.length-1, Math.ceil((ve-binOff)/binSz));
  if (firstBin > lastBin) return;
  const vd = data.slice(firstBin, lastBin+1);
  const maxCov = Math.max(1, ...vd), logMax = Math.log10(maxCov+1);

  const toH  = val => val<=0 ? 0 : Math.log10(val+1)/logMax * (trackH-4);
  const bx1  = i => (binOff+(firstBin+i)*binSz-vs)*scaleX;
  const bx2  = i => (binOff+(firstBin+i)*binSz+binSz-vs)*scaleX;
  const bcx  = i => (bx1(i)+bx2(i))/2;

  if (style === 'histogram') {
    for (let i=0; i<vd.length; i++) {
      if (!vd[i]) continue;
      const rx1=bx1(i), rx2=bx2(i); if (rx2<=rx1) continue;
      const fullW=rx2-rx1;
      const barW=Math.max(1, fullW*0.55)|0;
      const px=(rx1+(fullW-barW)*0.5)|0;
      const clampX=Math.max(0,px), clampW=Math.min(barW,W-clampX);
      if (clampW<=0) continue;
      const h=toH(vd[i])|0;
      const y = isPlus ? yBase+trackH-2-h : yBase+2;
      out.push(`<rect x="${clampX}" y="${y}" width="${clampW}" height="${h}" fill="${col}" opacity="0.75"/>`);
    }
  } else if (style === 'kernel') {
    const sm = gaussSmooth1D(vd, 0.8);
    const pts = isPlus
      ? sm.map((_,i)=>`${bcx(i)|0},${(yBase+trackH-2-toH(sm[i]))|0}`).join(' ')
      : sm.map((_,i)=>`${bcx(i)|0},${(yBase+2+toH(sm[i]))|0}`).join(' ');
    const base = isPlus ? yBase+trackH-2 : yBase+2;
    if (pts) {
      out.push(`<polygon points="${pts} ${W},${base} 0,${base}" fill="${col}" opacity="0.2"/>`);
      out.push(`<polyline points="${pts}" fill="none" stroke="${col}" stroke-width="1.5" opacity="0.9"/>`);
    }
  } else { // reads
    const actualReads = isPlus ? (state.reads_plus || []) : (state.reads_minus || []);
    if (actualReads.length > 0) {
      // Draw actual reads at exact genomic coordinates, stacked greedily
      const rh=2, rg=1, maxRows=Math.floor((trackH-6)/(rh+rg));
      const rowEnds=[];  // rightmost pixel of last read placed in each row
      for (const [rStart, rEnd] of actualReads) {
        const px1=Math.max(0, ((rStart-vs)*scaleX)|0);
        const px2=Math.min(W, Math.max(px1+1, ((rEnd-vs)*scaleX+0.5)|0));
        if (px2<=0 || px1>=W) continue;
        let row=0;
        while (row<rowEnds.length && rowEnds[row]>px1) row++;
        if (row>=maxRows) continue;
        rowEnds[row]=px2;
        const ry=isPlus ? yBase+trackH-4-row*(rh+rg) : yBase+3+row*(rh+rg);
        if (isPlus && ry<yBase+2) continue;
        if (!isPlus && ry+rh>yBase+trackH-2) continue;
        out.push(`<rect x="${px1}" y="${ry}" width="${px2-px1}" height="${rh}" fill="${col}" opacity="0.8" rx="0.5"/>`);
      }
    } else {
      // Fallback: simulate reads from bin coverage counts
      const rh=2, rg=1, maxR=Math.floor((trackH-6)/(rh+rg));
      function lcg(s) { return (Math.imul(s,1664525)+1013904223)|0; }
      for (let i=0; i<vd.length; i++) {
        if (!vd[i]) continue;
        const fx1=bx1(i), fx2=bx2(i);
        const px1=Math.max(0,fx1|0), px2=Math.min(W,(fx2+0.5)|0);
        if (px2<=px1) continue;
        const bw=px2-px1;
        const rw=Math.max(2, Math.min(18, (bw*0.38)|0));
        let seed=(firstBin+i)*(isPlus?12345:67890);
        for (let r=0; r<Math.min(maxR,vd[i]); r++) {
          seed=lcg(seed);
          const room=Math.max(0,bw-rw-2);
          const rx=px1+1+(room>0 ? (((seed>>8)&0xff)/256*room)|0 : 0);
          if (rx+rw>px2) continue;
          const ry=isPlus ? yBase+trackH-4-r*(rh+rg) : yBase+3+r*(rh+rg);
          if (isPlus && ry<yBase+2) break;
          if (!isPlus && ry+rh>yBase+trackH-2) break;
          out.push(`<rect x="${rx}" y="${ry}" width="${rw}" height="${rh}" fill="${col}" opacity="0.65" rx="0.5"/>`);
        }
      }
    }
  }

  out.push(`<text x="3" y="${isPlus?yBase+trackH-3:yBase+9}" font-size="7" fill="#484f58" font-family="monospace">${isPlus?'cov+':'cov-'}</text>`);
  for (const lv of [1,10,100,1000,10000]) {
    if (lv > maxCov*1.5) break;
    const h = toH(lv);
    const y = isPlus ? yBase+trackH-2-h : yBase+2+h;
    if (isPlus && y < yBase+4) break;
    if (!isPlus && y > yBase+trackH-4) break;
    out.push(`<line x1="0" y1="${y|0}" x2="5" y2="${y|0}" stroke="#2d3547" stroke-width="1"/>`);
    out.push(`<text x="7" y="${(y+3)|0}" font-size="7" fill="#3d4a5e" font-family="monospace">${lv<1000?lv:lv/1000+'k'}</text>`);
  }
}
