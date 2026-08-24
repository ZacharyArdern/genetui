// ── Pyodide init ──────────────────────────────────────────────────────────────
async function initPyodide() {
  const badge=document.getElementById('py-badge');
  badge.style.display=''; badge.textContent='py\u2026';
  setFigStatus('loading Pyodide\u2026');
  try {
    pyodide = await loadPyodide({indexURL:'https://cdn.jsdelivr.net/pyodide/v0.26.4/full/'});
    setFigStatus('installing packages\u2026');
    await pyodide.loadPackage(['micropip','matplotlib','numpy','biopython']);
    await pyodide.runPythonAsync('import micropip\nawait micropip.install("dna_features_viewer", keep_going=True)');
    badge.textContent='py ready'; badge.className='badge ready';
    setFigStatus('ready \u2014 renders when scrolling stops');
    pyReady=true;
    document.getElementById('btn-render').disabled=false;
    if (lastState) scheduleAutoRender();
    if (typeof window.updateFancyLayer === 'function') window.updateFancyLayer();
  } catch(err) {
    badge.textContent='py error'; badge.style.background='#da3633';
    setFigStatus('Pyodide error');
    document.getElementById('fig-output').innerHTML=`<pre class="fig-err">${esc(String(err))}</pre>`;
  }
}
(function(){
  const s=document.createElement('script');
  s.src='https://cdn.jsdelivr.net/pyodide/v0.26.4/full/pyodide.js';
  s.onload=()=>initPyodide();
  s.onerror=()=>{
    document.getElementById('py-badge').style.display='';
    document.getElementById('py-badge').textContent='offline';
    setFigStatus('CDN unavailable');
  };
  document.head.appendChild(s);
})();

// ── Auto-render debounce ──────────────────────────────────────────────────────
function scheduleAutoRender() {
  if (!pyReady) return;
  clearTimeout(debounceTimer);
  debounceTimer=setTimeout(()=>triggerRender(false), 700);
}

// ── Python code builder ───────────────────────────────────────────────────────
function buildPyCode(vs,ve,fmt) {
  const forPng = fmt==='png' || fmt===true;
  if (!lastState) return null;
  const width     = parseFloat(document.getElementById('fig-width').value);
  const fs        = parseInt(document.getElementById('fig-fs').value);
  const showLbls  = document.getElementById('fig-show-labels').checked;
  const hideLong  = document.getElementById('fig-hide-long-labels').checked;
  const ruler     = document.getElementById('fig-ruler').checked;
  const white     = document.getElementById('fig-white').checked;
  const sixframe  = document.getElementById('opt-sixframe').checked;
  const stopCod   = document.getElementById('opt-stopcodons').checked;
  const hlCodonEl = document.getElementById('opt-codons');
  const hlCodonsArr = hlCodonEl ? hlCodonEl.value.toUpperCase().split(/[,\s]+/).map(s=>s.trim()).filter(s=>s.length===3) : [];
  const pyHlCodons = hlCodonsArr.length ? `{'${hlCodonsArr.join("','")}'}` : 'set()';
  const seqLen    = Math.max(1,ve-vs);
  const feats     = lastState.features.filter(f=>f.start<=ve&&f.end>=vs);
  const bg        = white?"'white'":"'#0d1117'";
  const fg        = white?"'#111'":"'#c9d1d9'";
  const dpi       = forPng?600:96;
  const fmt_      = forPng?'png':(fmt==='pdf'?'pdf':'svg');

  const covStyle    = document.getElementById('cov-style-fig').value;
  const covOverlay  = document.getElementById('fig-cov-layout').value === 'overlay';
  const coordDp     = parseInt(document.getElementById('fig-coord-dp').value);
  const covLogScale = document.getElementById('fig-cov-log').checked;
  const covYMaxRaw  = parseFloat(document.getElementById('fig-cov-ymax').value);
  const covYMax     = isNaN(covYMaxRaw) ? 0 : covYMaxRaw;
  // Use _codonColors (kept in sync with fig-codon-lines + picker overrides)
  const figCodonCols = window._codonColors && Object.keys(window._codonColors).length
    ? window._codonColors : {};
  const pyFigCodonCols = JSON.stringify(figCodonCols);
  const covBinSz   = lastState.coverage_bin_size  || 1000;
  const covBinOff  = lastState.coverage_bin_start || 0;
  // coverage_plus/minus are Vec<Vec<u32>> (one inner array per BAM track).
  const rawTracksPlus  = lastState.coverage_plus  || [];
  const rawTracksMinus = lastState.coverage_minus || [];
  const isCovNested = rawTracksPlus.length > 0 && Array.isArray(rawTracksPlus[0]);
  const tracksPlus  = isCovNested ? rawTracksPlus  : (rawTracksPlus.length  ? [rawTracksPlus]  : []);
  const tracksMinus = isCovNested ? rawTracksMinus : (rawTracksMinus.length ? [rawTracksMinus] : []);
  const bamLabels   = lastState.bam_labels || tracksPlus.map((_,i) => `BAM ${i+1}`);
  const activeMask       = lastState.active_bam_mask || tracksPlus.map(() => true);
  const activeIndices    = bamLabels.map((_,i)=>i).filter(i => activeMask[i] !== false);
  const activeBamColors  = activeIndices.map(i => window.getBamEffectiveColor ? window.getBamEffectiveColor(i, bamLabels[i]||'') : _DEFAULT_PAL[i%6]);
  const pyBamColors      = JSON.stringify(activeBamColors);
  const extraGffLabels   = lastState.extra_gff_labels || [];
  const extraGffFeatures = lastState.extra_gff_features || [];
  const gffOverlay       = lastState.gff_overlay !== false;

  let pyCovTracks = '[]', pyCovBinSz = covBinSz, pyCovBinOff = 0;
  if (covStyle !== 'none' && tracksPlus.length) {
    const activeTP = tracksPlus.filter((_, i) => activeMask[i] !== false);
    const activeTM = tracksMinus.filter((_, i) => activeMask[i] !== false);
    const activeLabels = bamLabels.filter((_, i) => activeMask[i] !== false);
    if (activeTP.length) {
      const refTrack = activeTP[0];
      const fb = Math.max(0, Math.floor((vs - covBinOff) / covBinSz));
      const lb = Math.min(refTrack.length - 1, Math.ceil((ve - covBinOff) / covBinSz));
      if (fb <= lb) {
        pyCovBinOff = covBinOff + fb * covBinSz;
        pyCovTracks = JSON.stringify(activeTP.map((tp, i) => ({
          p: tp.slice(fb, lb+1),
          m: (activeTM[i] || []).slice(fb, lb+1),
          l: (activeLabels[i] || `BAM ${i+1}`).replace(/\.bam$/i,'').slice(0, 20),
        })));
      }
    }
  }

  // Per-read data: cap per track to avoid Pyodide MemoryError when tokenizing large literals.
  const MAX_READS_PER_TRACK = 600;
  const rawRP = lastState.reads_plus  || [];
  const rawRM = lastState.reads_minus || [];
  const isReadsNested = rawRP.length > 0 && Array.isArray(rawRP[0]);
  const readsTracksP = isReadsNested ? rawRP : (rawRP.length ? [rawRP] : []);
  const readsTracksM = isReadsNested ? rawRM : (rawRM.length ? [rawRM] : []);
  let pyReadsTracks = '[]';
  if (covStyle === 'reads' || covStyle === 'barplot') {
    const activeRP = readsTracksP.filter((_, i) => activeMask[i] !== false);
    const activeRM = readsTracksM.filter((_, i) => activeMask[i] !== false);
    const activeLabels = bamLabels.filter((_, i) => activeMask[i] !== false);
    if (activeRP.length) {
      const filt = r => r[1] >= vs && r[0] <= ve;
      pyReadsTracks = JSON.stringify(activeRP.map((rp, i) => ({
        p: rp.filter(filt).slice(0, MAX_READS_PER_TRACK),
        m: (activeRM[i] || []).filter(filt).slice(0, MAX_READS_PER_TRACK),
        l: (activeLabels[i] || `BAM ${i+1}`).replace(/\.bam$/i,'').slice(0, 20),
      })));
    }
  }
  const hasCov = covStyle !== 'none' && (pyCovTracks !== '[]' || pyReadsTracks !== '[]');

  // Contig boundary info for per-contig axis labels and separator lines
  const _cnames = lastState.contig_names || [];
  const _cstarts = lastState.contig_starts || [];
  const _cends = lastState.contig_ends || [];
  const pyContigs = JSON.stringify(_cnames.map((n,i)=>({n,s:_cstarts[i],e:_cends[i]})));
  // Contig separator positions visible in this viewport (relative coords)
  const visContigSeps = _cstarts.slice(1).filter(s => s > vs && s < ve).map(s => s - vs);

  const charWidthIn = (fs / 11) * parseFloat(document.getElementById('fig-overflow-thresh').value);
  function labelFitsInFeature(f) {
    const visibleBp = Math.min(f.end, ve) - Math.max(f.start, vs);
    const featWidthIn = visibleBp / seqLen * width;
    return featWidthIn >= f.name.length * charWidthIn;
  }

  const pyFeats = feats.map(f=>{
    const lbl = showLbls && !(hideLong && !labelFitsInFeature(f)) ? f.name.replace(/\\/g,'\\\\').replace(/'/g,"\\'") : '';
    const col = f.noncoding ? NC_COLOR : f.color;
    const s  = Math.max(0, f.start-vs);
    const ee = Math.min(seqLen, f.end-vs);
    const si = f.strand==='+'?1:-1;
    const fr = f.strand==='+'? (f.start-1)%3 : (f.end-1)%3;
    return `    (${s},${ee},${si},${lbl?`'${lbl}'`:'None'},'${col}',${f.noncoding?'True':'False'},${fr})`;
  }).join(',\n');

  const pyExtraFeats = JSON.stringify(extraGffFeatures.map((ef, gi) => {
    const visFeats = ef.filter(f => f.start <= ve && f.end >= vs);
    const colOverride = window._extraGffColors && window._extraGffColors[gi];
    return {
      label: (extraGffLabels[gi] || `GFF ${gi+2}`).replace(/\\/g,'\\\\').replace(/'/g,"\\'").slice(0,30),
      feats: visFeats.map(f => {
        const lbl2 = showLbls && !(hideLong && !labelFitsInFeature(f))
          ? f.name.replace(/\\/g,'\\\\').replace(/'/g,"\\'") : '';
        return [Math.max(0,f.start-vs), Math.min(seqLen,f.end-vs),
                f.strand==='+'?1:-1, lbl2||'',
                colOverride || f.color, f.noncoding?1:0, (f.strand==='+'?(f.start-1)%3:(f.end-1)%3)];
      })
    };
  }));
  const pyGffOverlay = gffOverlay ? 'True' : 'False';

  if (sixframe) {
    return `
import io,base64,traceback,matplotlib,matplotlib.gridspec as gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
_result=None
try:
    from dna_features_viewer import GraphicRecord,GraphicFeature
    _bg=${bg}; _fg=${fg}
    _seqLen=${seqLen}; _width=${width}; _ruler=${ruler?'True':'False'}
    _stopCod=${stopCod?'True':'False'}
    _feats_raw=[
${pyFeats}
    ]
    _STOP={'TAA','TAG','TGA'}
    _COMP={'A':'T','T':'A','G':'C','C':'G'}
    def _rc(s):
        return ''.join(_COMP.get(c,'N') for c in reversed(s.upper()))
    _seq_sub=globals().get('_seq_sub',''); _seq_sub_off=globals().get('_seq_sub_off',0)
    _geno_off=${vs}+_seq_sub_off
    _fwd_stops=[[],[],[]]; _rev_stops=[[],[],[]]
    _HL_CODONS=${pyHlCodons}
    _fwd_hl=[[],[],[]]; _rev_hl=[[],[],[]]
    if _seq_sub:
        _u=_seq_sub.upper(); _rc_u=_rc(_seq_sub)
        for _i in range(len(_u)-2):
            _cdn=_u[_i:_i+3]; _gfr=(_geno_off+_i)%3
            if _stopCod and _cdn in _STOP: _fwd_stops[_gfr].append(_seq_sub_off+_i)
            if _HL_CODONS and _cdn in _HL_CODONS: _fwd_hl[_gfr].append(_seq_sub_off+_i)
        for _i in range(len(_rc_u)-2):
            _cdn=_rc_u[_i:_i+3]; _gend=_geno_off+len(_seq_sub)-_i-1; _gfr=_gend%3
            if _stopCod and _cdn in _STOP: _rev_stops[_gfr].append(_seq_sub_off+len(_seq_sub)-_i-2)
            if _HL_CODONS and _cdn in _HL_CODONS: _rev_hl[_gfr].append(_seq_sub_off+len(_seq_sub)-_i-2)
    _ROWS=[('+1',1,0),('+2',1,1),('+3',1,2),('-1',-1,0),('-2',-1,1),('-3',-1,2)]
    _cov_style='${covStyle}'; _cov_overlay=${covOverlay?'True':'False'}
    import json as _jcc
    _cov_tracks=_jcc.loads('${pyCovTracks.replace(/\\/g,'\\\\').replace(/'/g,"\\'")}'); _reads_tracks=_jcc.loads('${pyReadsTracks.replace(/\\/g,'\\\\').replace(/'/g,"\\'")}')
    _cov_bin_sz=${pyCovBinSz}; _cov_bin_off=${pyCovBinOff}
    _cov_log=${covLogScale?'True':'False'}; _cov_ymax=${covYMax}
    _FIG_CODON_COLS=_jcc.loads('${pyFigCodonCols.replace(/\\/g,'\\\\').replace(/'/g,"\\'")}' or '{}')
    _has_cov=_cov_style!='none' and (len(_cov_tracks)>0 or len(_reads_tracks)>0)
    _extra_gffs=_jcc.loads('${pyExtraFeats.replace(/\\/g,'\\\\').replace(/'/g,"\\'")}'); _gff_overlay=${pyGffOverlay}
    if _extra_gffs and _gff_overlay:
        for _eg in _extra_gffs:
            for _ef in _eg['feats']:
                _feats_raw.append((_ef[0],_ef[1],_ef[2],_ef[3],_ef[4],bool(_ef[5]),_ef[6]))
    _TRACK_COLS_P=${pyBamColors} or ['#3fc9d0','#56d364','#a371f7','#f85149','#e3b341','#d66adc']
    _TRACK_COLS_M=_TRACK_COLS_P
    # row order: cov+? / nc+ / +1/+2/+3 / sep / -1/-2/-3 / nc- / cov-?
    _hr=([2.8] if _has_cov else [])+[0.32,0.52,0.52,0.52,0.25,0.52,0.52,0.52,0.32]+([2.8] if _has_cov else [])
    _off=1 if _has_cov else 0
    _fig_h=6*0.20+2*0.15+0.75+(3.2 if _has_cov else 0)
    fig=plt.figure(figsize=(_width,_fig_h))
    fig.patch.set_facecolor(_bg)
    _nt_bam=len(_cov_tracks) if _cov_style not in ('reads','barplot') else len(_reads_tracks)
    _fig_right=0.75 if _nt_bam>1 else 0.98
    gs=gridspec.GridSpec(len(_hr),1,figure=fig,hspace=0.0,
        height_ratios=_hr,top=0.96,bottom=0.09,left=0.11,right=_fig_right)
    import matplotlib.ticker as _tck,json as _json
    _CONTIGS=_json.loads('${pyContigs}')
    _CONTIG_SEPS=${JSON.stringify(visContigSeps)}
    _coord_dp=${coordDp}; _span=${seqLen}
    def _mb_dp(): return _coord_dp if _coord_dp else (4 if _span<10000 else 3 if _span<100000 else 2)
    def _kb_dp(): return _coord_dp if _coord_dp else (3 if _span<10000 else 2 if _span<100000 else 1)
    def _gfmt(x,p):
        g=int(x)+${vs}
        if _CONTIGS:
            for _c in _CONTIGS:
                if _c['s']<=g<=_c['e']:
                    _v=g-_c['s']+1
                    if _v>=1000000: return f'{_v/1000000:.{_mb_dp()}f}M'
                    if _v>=1000: return f'{_v/1000:.{_kb_dp()}f}k'
                    return str(_v)
        _v=g
        if _v>=1000000: return f'{_v/1000000:.{_mb_dp()}f}M'
        if _v>=1000: return f'{_v/1000:.{_kb_dp()}f}k'
        return str(_v)
    _bpfmt=_tck.FuncFormatter(_gfmt)
    _axes=[]
    for _ri,(_lbl,_si,_fr) in enumerate(_ROWS):
        # nc+ occupies _off+0, so +frames start at _off+1; separator at _off+4; -frames at _off+5..7
        _gsi = _off+1+_ri if _ri<3 else _off+_ri+2
        ax=fig.add_subplot(gs[_gsi])
        _gfs=[]
        for (s,e,st,lb,col,noncoding,frame) in _feats_raw:
            if st!=_si or noncoding: continue
            if frame!=_fr: continue
            _gfs.append(GraphicFeature(start=s,end=e,strand=_si,label=(lb if lb else None),color=col))
        rec=GraphicRecord(sequence_length=_seqLen,features=_gfs)
        rec.plot(ax=ax,with_ruler=False,draw_line=True)
        ax.set_ylim(-0.55,0.55)
        ax.set_facecolor(_bg)
        ax.set_ylabel(_lbl,fontsize=7,rotation=0,labelpad=22,va='center',color=_fg)
        for sp in ax.spines.values(): sp.set_color(_fg)
        ax.tick_params(colors=_fg,labelsize=${fs-2})
        ax.set_xticklabels([])
        if _seq_sub:
            if _stopCod:
                for _x in (_fwd_stops[_fr] if _si==1 else _rev_stops[_fr]):
                    ax.axvline(_x,color='black',alpha=0.7,linewidth=0.8)
            if _HL_CODONS:
                for _x in (_fwd_hl[_fr] if _si==1 else _rev_hl[_fr]):
                    ax.axvspan(_x,_x+3,color='#3fb950',alpha=0.32,linewidth=0,zorder=3)
            for _sx in _CONTIG_SEPS:
                ax.axvline(_sx,color='#58a6ff',linewidth=0.8,linestyle='--',alpha=0.6,zorder=4)
        _axes.append(ax)
    _nc_p=[GraphicFeature(start=s,end=e,strand=st,label=(lb if lb else None),color=col) for (s,e,st,lb,col,nc,fr) in _feats_raw if nc and st==1]
    _nc_m=[GraphicFeature(start=s,end=e,strand=st,label=(lb if lb else None),color=col) for (s,e,st,lb,col,nc,fr) in _feats_raw if nc and st==-1]
    for (_nc_feats,_nc_lbl,_nc_idx,_show_x) in [(_nc_p,'nc+',_off+0,False),(_nc_m,'nc-',_off+8,not _has_cov)]:
        _anc=fig.add_subplot(gs[_nc_idx])
        GraphicRecord(sequence_length=_seqLen,features=_nc_feats).plot(ax=_anc,with_ruler=_show_x,draw_line=True)
        _anc.set_ylim(-0.55,0.55)
        _anc.set_facecolor(_bg)
        _anc.set_ylabel(_nc_lbl,fontsize=7,rotation=0,labelpad=22,va='center',color=_fg)
        for sp in _anc.spines.values(): sp.set_color(_fg)
        _anc.tick_params(colors=_fg,labelsize=${fs-2})
        if _show_x:
            _anc.xaxis.set_major_formatter(_bpfmt)
            _anc.tick_params(axis='x',colors=_fg,labelsize=${fs-2})
        for _sx in _CONTIG_SEPS:
            _anc.axvline(_sx,color='#58a6ff',linewidth=0.8,linestyle='--',alpha=0.6,zorder=4)
    if _FIG_CODON_COLS and _seq_sub:
        _u=_seq_sub.upper(); _rc_u2=_rc(_seq_sub)
        _fcd_fwd=[{} for _ in range(3)]; _fcd_rev=[{} for _ in range(3)]
        for _i in range(len(_u)-2):
            _cdn=_u[_i:_i+3]
            if _cdn in _FIG_CODON_COLS:
                _fcd_fwd[(_geno_off+_i)%3].setdefault(_cdn,[]).append(_seq_sub_off+_i)
        for _i in range(len(_rc_u2)-2):
            _cdn=_rc_u2[_i:_i+3]
            if _cdn in _FIG_CODON_COLS:
                _fcd_rev[(_geno_off+len(_seq_sub)-_i-1)%3].setdefault(_cdn,[]).append(_seq_sub_off+len(_seq_sub)-_i-2)
        for _ri2,(_lbl2,_si2,_fr2) in enumerate(_ROWS):
            _fd=_fcd_fwd[_fr2] if _si2==1 else _fcd_rev[_fr2]
            for _cdn,_pxs in _fd.items():
                _cc=_FIG_CODON_COLS[_cdn]
                for _x in _pxs: _axes[_ri2].axvline(_x,color=_cc,linewidth=0.8,alpha=0.7,zorder=3)
    _ax_sep=fig.add_subplot(gs[_off+4])
    _ax_sep.set_visible(False)
    _pos3 =_axes[2].get_position()
    _pos4 =_axes[3].get_position()
    _ymid =(_pos3.y0+_pos4.y1)/2
    from matplotlib.lines import Line2D
    fig.add_artist(Line2D([0.07,0.98],[_ymid,_ymid],transform=fig.transFigure,
        color=_fg,linewidth=0.8,linestyle='-'))
    if _has_cov:
        import numpy as _np
        from matplotlib.patches import Rectangle; from matplotlib.collections import PatchCollection
        _bin_len=len(_cov_tracks[0]['p']) if _cov_tracks else 0
        _xs=[(_cov_bin_off+i*_cov_bin_sz-${vs}) for i in range(_bin_len)]
        def _cov_panel(ax,strand,invert=False,show_xaxis=False):
            _tcols=_TRACK_COLS_P if strand=='+' else _TRACK_COLS_M
            _tdata=_cov_tracks if _cov_style not in ('reads','barplot') else []
            _rdata=_reads_tracks if _cov_style in ('reads','barplot') else []
            _nt=max(len(_tdata),len(_rdata),1)
            ax.set_facecolor(_bg); ax.set_xlim(0,${seqLen})
            for sp in ax.spines.values(): sp.set_color(_fg)
            ax.tick_params(colors=_fg,labelsize=${fs-2})
            _src=_tdata or _rdata
            _lbl_src=_cov_tracks or _reads_tracks
            _ylbl=_lbl_src[0].get('l','')[:16] if len(_lbl_src)==1 else ('cov+' if strand=='+' else 'cov-')
            ax.set_ylabel(_ylbl,fontsize=6,rotation=0,labelpad=38,va='center',color=_fg)
            _rh,_rg,_mr=0.8,0.3,25
            if _cov_style=='reads':
                _min_rw=max(float(${seqLen})/(${width}*${dpi})*2,1.0)
                for _ti,_td in enumerate(_rdata):
                    _col=_tcols[_ti%len(_tcols)]
                    _reads_list=_td['p'] if strand=='+' else _td['m']
                    _rects=[]; _row_ends=[]
                    _row_off=0 if _cov_overlay else ((_ti*(_mr//_nt)) if _nt>1 else 0)
                    _row_cap=_mr if _cov_overlay else _mr//_nt
                    for (_rs,_re) in sorted(_reads_list,key=lambda r:r[0]):
                        _rx=_rs-${vs}; _rw=max(_min_rw,float(_re-_rs))
                        if _rx+_rw<0 or _rx>${seqLen}: continue
                        _row=0
                        while _row<len(_row_ends) and _row_ends[_row]>_rx: _row+=1
                        if _row>=_row_cap: continue
                        if _row>=len(_row_ends): _row_ends.append(_rx+_rw)
                        else: _row_ends[_row]=_rx+_rw
                        _rects.append(Rectangle((_rx,(_row+_row_off)*(_rh+_rg)),_rw,_rh,label=_td['l']))
                    _alp=0.45 if (_cov_overlay and _nt>1) else 0.8
                    if _rects: ax.add_collection(PatchCollection(_rects,facecolor=_col,alpha=_alp,linewidth=0))
                ax.set_ylim(0,_mr*(_rh+_rg)); ax.set_xlim(0,${seqLen})
                ax.set_yticks([]); ax.set_yticklabels([])
                if _nt>1:
                    from matplotlib.patches import Patch
                    _hnd=[Patch(facecolor=_tcols[i%len(_tcols)],label=_rdata[i]['l']) for i in range(len(_rdata))]
                    ax.legend(handles=_hnd,fontsize=5,bbox_to_anchor=(1.02,1),loc='upper left',borderaxespad=0,framealpha=0.5,facecolor=_bg,edgecolor=_fg,labelcolor=_fg)
            elif _cov_style=='barplot':
                _min_rw=max(float(${seqLen})/(${width}*${dpi})*2,1.0)
                _max_h=0
                for _ti,_td in enumerate(_rdata):
                    _col=_tcols[_ti%len(_tcols)]
                    _reads_list=_td['p'] if strand=='+' else _td['m']
                    _grps={}
                    for (_rs,_re) in _reads_list:
                        _rx=float(_rs-${vs}); _rw=float(_re-_rs)
                        if _rx+_rw<0 or _rx>${seqLen}: continue
                        _rx=max(0.0,_rx)
                        _k=(round(_rx),max(1,round(_rw)))
                        _grps[_k]=_grps.get(_k,0)+1
                    _by_x={}
                    for (_xk,_wk),_cnt in _grps.items():
                        _by_x.setdefault(_xk,[]).append((_wk,_cnt))
                    _ti_rects=[]
                    for _xk,_stk in _by_x.items():
                        _y=0
                        for _wk,_cnt in sorted(_stk,reverse=True):
                            _ti_rects.append(Rectangle((_xk,_y),max(float(_wk),_min_rw),float(_cnt)))
                            _y+=_cnt
                        _max_h=max(_max_h,_y)
                    _alp=0.45 if (_cov_overlay and _nt>1) else 0.8
                    if _ti_rects: ax.add_collection(PatchCollection(_ti_rects,facecolor=_col,alpha=_alp,linewidth=0))
                ax.set_ylim(0,max(_max_h,1)*1.15); ax.set_xlim(0,${seqLen})
                ax.tick_params(axis='y',colors=_fg,labelsize=6)
                if _nt>1:
                    from matplotlib.patches import Patch
                    _hnd=[Patch(facecolor=_tcols[i%len(_tcols)],label=_rdata[i]['l']) for i in range(len(_rdata))]
                    ax.legend(handles=_hnd,fontsize=5,bbox_to_anchor=(1.02,1),loc='upper left',borderaxespad=0,framealpha=0.5,facecolor=_bg,edgecolor=_fg,labelcolor=_fg)
            else:
                _bw=_cov_bin_sz*0.9/(1 if _cov_overlay else max(_nt,1)); _all_mx=[]
                for _ti,_td in enumerate(_tdata):
                    _col=_tcols[_ti%len(_tcols)]
                    _raw=_td['p'] if strand=='+' else _td['m']
                    lv=_np.log10(_np.array(_raw,dtype=float)+1)
                    _all_mx.append(float(lv.max()) if len(lv) else 0)
                    _xoff=0 if _cov_overlay else (_ti*_bw if _nt>1 else 0)
                    _alp=0.5 if (_cov_overlay and _nt>1) else 0.75
                    if _cov_style=='histogram':
                        ax.bar([x+_xoff for x in _xs],list(lv),width=_bw,align='edge',color=_col,alpha=_alp,label=_td['l'])
                    else:
                        _sig=0.8; _rr=int(_np.ceil(3*_sig))
                        _k=_np.exp(-0.5*(_np.arange(-_rr,_rr+1)/_sig)**2); _k/=_k.sum()
                        _sm=_np.convolve(lv,_k,mode='same')
                        ax.fill_between(_xs,_sm,alpha=(0.12 if _nt>1 else 0.2),color=_col)
                        ax.plot(_xs,_sm,color=_col,linewidth=1.2,label=_td['l'])
                _mx=max(_all_mx) if _all_mx else 0.01
                ax.set_ylim(0,max(_mx,0.01)*1.1)
                _ytv=[v for v in [1,10,100,1000,10000] if _np.log10(v+1)<=_mx*1.05]
                ax.set_yticks([_np.log10(v+1) for v in _ytv]); ax.set_yticklabels([str(v) for v in _ytv],fontsize=6)
                if _nt>1: ax.legend(fontsize=5,bbox_to_anchor=(1.02,1),loc='upper left',borderaxespad=0,framealpha=0.5,facecolor=_bg,edgecolor=_fg,labelcolor=_fg)
            if invert: ax.invert_yaxis()
            if show_xaxis:
                ax.xaxis.set_major_formatter(_bpfmt); ax.tick_params(axis='x',colors=_fg,labelsize=${fs-2})
            else:
                ax.set_xticklabels([])
        _cov_panel(fig.add_subplot(gs[0]),'+')
        _cov_panel(fig.add_subplot(gs[-1]),'-',invert=True,show_xaxis=True)
    buf=io.BytesIO()
    fig.savefig(buf,format='${fmt_}',dpi=${dpi},bbox_inches='tight',facecolor=_bg)
    plt.close(fig); buf.seek(0)
    _result='OK:'+base64.b64encode(buf.read()).decode()
except Exception:
    _result='ERR:'+traceback.format_exc()
`;
  }

  // Simple mode — features as tuples (s,e,si,lbl,col,noncoding)
  const pyFeatsSingle = feats.map(f=>{
    const lbl=showLbls&&!(hideLong&&!labelFitsInFeature(f))?f.name.replace(/\\/g,'\\\\').replace(/'/g,"\\'"):'';
    const col=f.noncoding?NC_COLOR:f.color;
    const s=Math.max(0,f.start-vs), ee=Math.min(seqLen,f.end-vs);
    const si=f.strand==='+'?1:-1;
    const isp=f.strand==='+'?'True':'False';
    return `    (${s},${ee},${si},${lbl?`'${lbl}'`:'None'},'${col}',${f.noncoding?'True':'False'},${isp})`;
  }).join(',\n');
  return `
import io,base64,traceback,matplotlib,matplotlib.gridspec as gridspec,matplotlib.ticker as _tck
matplotlib.use('Agg')
import matplotlib.pyplot as plt
_result=None
try:
    from dna_features_viewer import GraphicRecord,GraphicFeature
    _bg=${bg}; _fg=${fg}; _seqLen=${seqLen}; _width=${width}
    _feats_raw=[
${pyFeatsSingle}
    ]
    import json as _json
    _CONTIGS=_json.loads('${pyContigs}')
    _CONTIG_SEPS=${JSON.stringify(visContigSeps)}
    _coord_dp=${coordDp}; _span=${seqLen}
    def _mb_dp(): return _coord_dp if _coord_dp else (4 if _span<10000 else 3 if _span<100000 else 2)
    def _kb_dp(): return _coord_dp if _coord_dp else (3 if _span<10000 else 2 if _span<100000 else 1)
    def _gfmt(x,p):
        g=int(x)+${vs}
        if _CONTIGS:
            for _c in _CONTIGS:
                if _c['s']<=g<=_c['e']:
                    _v=g-_c['s']+1
                    if _v>=1000000: return f'{_v/1000000:.{_mb_dp()}f}M'
                    if _v>=1000: return f'{_v/1000:.{_kb_dp()}f}k'
                    return str(_v)
        _v=g
        if _v>=1000000: return f'{_v/1000000:.{_mb_dp()}f}M'
        if _v>=1000: return f'{_v/1000:.{_kb_dp()}f}k'
        return str(_v)
    _bpfmt=_tck.FuncFormatter(_gfmt)
    _COMP={'A':'T','T':'A','G':'C','C':'G'}
    def _rc(s): return ''.join(_COMP.get(c,'N') for c in reversed(s.upper()))
    _seq_sub=globals().get('_seq_sub',''); _seq_sub_off=globals().get('_seq_sub_off',0)
    _geno_off=${vs}+_seq_sub_off
    _HL_CODONS=${pyHlCodons}
    _fwd_hl=[]; _rev_hl=[]
    if _HL_CODONS and _seq_sub:
        _u=_seq_sub.upper(); _rc_u=_rc(_seq_sub)
        for _i in range(len(_u)-2):
            if _u[_i:_i+3] in _HL_CODONS: _fwd_hl.append(_seq_sub_off+_i)
        for _i in range(len(_rc_u)-2):
            if _rc_u[_i:_i+3] in _HL_CODONS: _rev_hl.append(_seq_sub_off+len(_seq_sub)-_i-2)
    _cov_style='${covStyle}'; _cov_overlay=${covOverlay?'True':'False'}
    import json as _jcc
    _cov_tracks=_jcc.loads('${pyCovTracks.replace(/\\/g,'\\\\').replace(/'/g,"\\'")}'); _reads_tracks=_jcc.loads('${pyReadsTracks.replace(/\\/g,'\\\\').replace(/'/g,"\\'")}')
    _cov_bin_sz=${pyCovBinSz}; _cov_bin_off=${pyCovBinOff}
    _cov_log=${covLogScale?'True':'False'}; _cov_ymax=${covYMax}
    _FIG_CODON_COLS=_jcc.loads('${pyFigCodonCols.replace(/\\/g,'\\\\').replace(/'/g,"\\'")}' or '{}')
    _has_cov=_cov_style!='none' and (len(_cov_tracks)>0 or len(_reads_tracks)>0)
    _extra_gffs=_jcc.loads('${pyExtraFeats.replace(/\\/g,'\\\\').replace(/'/g,"\\'")}'); _gff_overlay=${pyGffOverlay}
    if _extra_gffs and _gff_overlay:
        for _eg in _extra_gffs:
            for _ef in _eg['feats']:
                _feats_raw.append((_ef[0],_ef[1],_ef[2],_ef[3],_ef[4],bool(_ef[5]),bool(_ef[2]==1)))
    _f_c   =[GraphicFeature(start=s,end=e,strand=si,label=(lb if lb else None),color=col) for (s,e,si,lb,col,nc,isp) in _feats_raw if not nc]
    _f_nc_p=[GraphicFeature(start=s,end=e,strand=si,label=(lb if lb else None),color=col) for (s,e,si,lb,col,nc,isp) in _feats_raw if nc and isp]
    _f_nc_m=[GraphicFeature(start=s,end=e,strand=si,label=(lb if lb else None),color=col) for (s,e,si,lb,col,nc,isp) in _feats_raw if nc and not isp]
    _has_c=bool(_f_c); _has_nc_p=bool(_f_nc_p); _has_nc_m=bool(_f_nc_m)
    _TRACK_COLS_P=${pyBamColors} or ['#3fc9d0','#56d364','#a371f7','#f85149','#e3b341','#d66adc']
    _TRACK_COLS_M=_TRACK_COLS_P
    # Assign explicit gridspec row indices: cov+ / nc+ / coding / nc- / cov- / extra GFF tracks
    _nrows=0; _hr=[]
    _ri_cov_p=_ri_nc_p=_ri_c=_ri_nc_m=_ri_cov_m=-1
    if _has_cov:  _ri_cov_p=_nrows; _nrows+=1; _hr.append(2.8)
    if _has_nc_p: _ri_nc_p =_nrows; _nrows+=1; _hr.append(0.46)
    if _has_c:    _ri_c    =_nrows; _nrows+=1; _hr.append(1.75)
    if _has_nc_m: _ri_nc_m =_nrows; _nrows+=1; _hr.append(0.46)
    if _has_cov:  _ri_cov_m=_nrows; _nrows+=1; _hr.append(2.8)
    _extra_gff_rows=[]
    if _extra_gffs and not _gff_overlay:
        for _eg in _extra_gffs:
            _extra_gff_rows.append(_nrows); _nrows+=1; _hr.append(0.63)
    _fig_h=sum(_hr)+0.55
    _nt_bam=len(_cov_tracks) if _cov_style not in ('reads','barplot') else len(_reads_tracks)
    _fig_right=0.75 if _nt_bam>1 else 0.98
    fig=plt.figure(figsize=(_width,_fig_h)); fig.patch.set_facecolor(_bg)
    gs=gridspec.GridSpec(_nrows,1,figure=fig,height_ratios=_hr,hspace=0.04,top=0.97,bottom=0.09,left=0.11,right=_fig_right)
    import numpy as _np
    from matplotlib.patches import Rectangle; from matplotlib.collections import PatchCollection
    _bin_len=len(_cov_tracks[0]['p']) if _cov_tracks else 0
    _xs=[(_cov_bin_off+i*_cov_bin_sz-${vs}) for i in range(_bin_len)] if _has_cov else []
    def _cov_panel(ax,strand,invert=False,show_xaxis=False):
        _tcols=_TRACK_COLS_P if strand=='+' else _TRACK_COLS_M
        _tdata=_cov_tracks if _cov_style not in ('reads','barplot') else []
        _rdata=_reads_tracks if _cov_style in ('reads','barplot') else []
        _nt=max(len(_tdata),len(_rdata),1)
        ax.set_facecolor(_bg); ax.set_xlim(0,_seqLen)
        for sp in ax.spines.values(): sp.set_color(_fg)
        ax.tick_params(colors=_fg,labelsize=${fs-2})
        _src=_tdata or _rdata
        _lbl_src=_cov_tracks or _reads_tracks
        _ylbl=_lbl_src[0].get('l','')[:16] if len(_lbl_src)==1 else ('cov+' if strand=='+' else 'cov-')
        ax.set_ylabel(_ylbl,fontsize=6,rotation=0,labelpad=38,va='center',color=_fg)
        _rh,_rg,_mr=0.8,0.3,25
        if _cov_style=='reads':
            _min_rw=max(_seqLen/(${width}*${dpi})*2,1.0)
            for _ti,_td in enumerate(_rdata):
                _col=_tcols[_ti%len(_tcols)]
                _reads_list=_td['p'] if strand=='+' else _td['m']
                _rects=[]; _row_ends=[]
                _row_off=0 if _cov_overlay else ((_ti*(_mr//_nt)) if _nt>1 else 0)
                _row_cap=_mr if _cov_overlay else _mr//_nt
                for (_rs,_re) in sorted(_reads_list,key=lambda r:r[0]):
                    _rx=_rs-${vs}; _rw=max(_min_rw,float(_re-_rs))
                    if _rx+_rw<0 or _rx>_seqLen: continue
                    _row=0
                    while _row<len(_row_ends) and _row_ends[_row]>_rx: _row+=1
                    if _row>=_row_cap: continue
                    if _row>=len(_row_ends): _row_ends.append(_rx+_rw)
                    else: _row_ends[_row]=_rx+_rw
                    _rects.append(Rectangle((_rx,(_row+_row_off)*(_rh+_rg)),_rw,_rh))
                _alp=0.45 if (_cov_overlay and _nt>1) else 0.8
                if _rects: ax.add_collection(PatchCollection(_rects,facecolor=_col,alpha=_alp,linewidth=0))
            ax.set_ylim(0,_mr*(_rh+_rg)); ax.set_xlim(0,_seqLen)
            ax.set_yticks([]); ax.set_yticklabels([])
            if _nt>1:
                from matplotlib.patches import Patch
                _hnd=[Patch(facecolor=_tcols[i%len(_tcols)],label=_rdata[i]['l']) for i in range(len(_rdata))]
                ax.legend(handles=_hnd,fontsize=5,bbox_to_anchor=(1.02,1),loc='upper left',borderaxespad=0,framealpha=0.5,facecolor=_bg,edgecolor=_fg,labelcolor=_fg)
        elif _cov_style=='barplot':
            _min_rw=max(float(_seqLen)/(${width}*${dpi})*2,1.0)
            _max_h=0
            for _ti,_td in enumerate(_rdata):
                _col=_tcols[_ti%len(_tcols)]
                _reads_list=_td['p'] if strand=='+' else _td['m']
                _grps={}
                for (_rs,_re) in _reads_list:
                    _rx=float(_rs-${vs}); _rw=float(_re-_rs)
                    if _rx+_rw<0 or _rx>_seqLen: continue
                    _rx=max(0.0,_rx)
                    _k=(round(_rx),max(1,round(_rw)))
                    _grps[_k]=_grps.get(_k,0)+1
                _by_x={}
                for (_xk,_wk),_cnt in _grps.items():
                    _by_x.setdefault(_xk,[]).append((_wk,_cnt))
                _ti_rects=[]
                for _xk,_stk in _by_x.items():
                    _y=0
                    for _wk,_cnt in sorted(_stk,reverse=True):
                        _ti_rects.append(Rectangle((_xk,_y),max(float(_wk),_min_rw),float(_cnt)))
                        _y+=_cnt
                    _max_h=max(_max_h,_y)
                _alp=0.45 if (_cov_overlay and _nt>1) else 0.8
                if _ti_rects: ax.add_collection(PatchCollection(_ti_rects,facecolor=_col,alpha=_alp,linewidth=0))
            ax.set_ylim(0,max(_max_h,1)*1.15); ax.set_xlim(0,_seqLen)
            ax.tick_params(axis='y',colors=_fg,labelsize=6)
            if _nt>1:
                from matplotlib.patches import Patch
                _hnd=[Patch(facecolor=_tcols[i%len(_tcols)],label=_rdata[i]['l']) for i in range(len(_rdata))]
                ax.legend(handles=_hnd,fontsize=5,bbox_to_anchor=(1.02,1),loc='upper left',borderaxespad=0,framealpha=0.5,facecolor=_bg,edgecolor=_fg,labelcolor=_fg)
        else:
            _bw=_cov_bin_sz*0.9/(1 if _cov_overlay else max(_nt,1)); _all_mx=[]
            for _ti,_td in enumerate(_tdata):
                _col=_tcols[_ti%len(_tcols)]
                _raw=_td['p'] if strand=='+' else _td['m']
                _ra=_np.array(_raw,dtype=float); lv=(_np.log10(_ra+1) if _cov_log else _ra)
                _all_mx.append(float(lv.max()) if len(lv) else 0)
                _xoff=0 if _cov_overlay else (_ti*_bw if _nt>1 else 0)
                _alp=0.5 if (_cov_overlay and _nt>1) else 0.75
                if _cov_style=='histogram':
                    ax.bar([x+_xoff for x in _xs],list(lv),width=_bw,align='edge',color=_col,alpha=_alp,label=_td['l'])
                else:
                    _sig=0.8; _rr=int(_np.ceil(3*_sig))
                    _k=_np.exp(-0.5*(_np.arange(-_rr,_rr+1)/_sig)**2); _k/=_k.sum()
                    _sm=_np.convolve(lv,_k,mode='same')
                    ax.fill_between(_xs,_sm,alpha=(0.12 if _nt>1 else 0.2),color=_col)
                    ax.plot(_xs,_sm,color=_col,linewidth=1.2,label=_td['l'])
            _mx=max(_all_mx) if _all_mx else 0.01
            _ym=_cov_ymax if _cov_ymax>0 else max(_mx,0.01)*1.1
            ax.set_ylim(0,_ym)
            if _cov_log:
                _ytv=[v for v in [1,10,100,1000,10000] if _np.log10(v+1)<=_ym*1.05]
                ax.set_yticks([_np.log10(v+1) for v in _ytv]); ax.set_yticklabels([str(v) for v in _ytv],fontsize=6)
            else:
                ax.tick_params(axis='y',colors=_fg,labelsize=6)
            if _nt>1: ax.legend(fontsize=5,bbox_to_anchor=(1.02,1),loc='upper left',borderaxespad=0,framealpha=0.5,facecolor=_bg,edgecolor=_fg,labelcolor=_fg)
        if invert: ax.invert_yaxis()
        if show_xaxis:
            ax.xaxis.set_major_locator(_tck.AutoLocator()); ax.xaxis.set_major_formatter(_bpfmt)
            ax.tick_params(axis='x',colors=_fg,labelsize=${fs-2})
        else:
            ax.set_xticklabels([])
    def _gene_panel(ax,rec,label,lpad,show_xaxis):
        rec.plot(ax=ax,with_ruler=show_xaxis,draw_line=True)
        ax.set_facecolor(_bg)
        for sp in ax.spines.values(): sp.set_color(_fg)
        ax.set_ylabel(label,fontsize=7,rotation=0,labelpad=lpad,va='center',color=_fg)
        for t in ax.texts: t.set_fontsize(${fs})
        if show_xaxis:
            ax.xaxis.set_major_formatter(_bpfmt)
            ax.tick_params(axis='x',colors=_fg,labelsize=${fs-2})
        return ax
    _has_extra_tracks = bool(_extra_gff_rows)
    _nc_p_show_x =_has_nc_p and not _has_c and not _has_nc_m and not _has_cov and not _has_extra_tracks
    _c_show_x    =_has_c    and not _has_nc_m and not _has_cov and not _has_extra_tracks
    _nc_m_show_x =_has_nc_m and not _has_cov and not _has_extra_tracks
    _g_axes=[]
    if _ri_cov_p>=0: _cov_panel(fig.add_subplot(gs[_ri_cov_p]),'+')
    if _ri_nc_p >=0: _g_axes.append(_gene_panel(fig.add_subplot(gs[_ri_nc_p]), GraphicRecord(sequence_length=_seqLen,features=_f_nc_p),'nc+',22,_nc_p_show_x))
    if _ri_c    >=0: _g_axes.append(_gene_panel(fig.add_subplot(gs[_ri_c]),    GraphicRecord(sequence_length=_seqLen,features=_f_c),    'genes',30,_c_show_x))
    if _ri_nc_m >=0: _g_axes.append(_gene_panel(fig.add_subplot(gs[_ri_nc_m]), GraphicRecord(sequence_length=_seqLen,features=_f_nc_m),'nc-',22,_nc_m_show_x))
    if _ri_cov_m>=0: _cov_panel(fig.add_subplot(gs[_ri_cov_m]),'-',invert=True,show_xaxis=(not _has_extra_tracks))
    for _egi,_eri in enumerate(_extra_gff_rows):
        _eg=_extra_gffs[_egi]
        _egfs=[GraphicFeature(start=max(0,s),end=min(_seqLen,max(s+1,e)),strand=si,label=(lb if lb else None),color=col) for (s,e,si,lb,col,nc,fr) in _eg['feats'] if e>s]
        _is_last=(_egi==len(_extra_gff_rows)-1)
        if not _egfs:
            _ax_eg=fig.add_subplot(gs[_eri]); _ax_eg.set_visible(False)
        else:
            try:
                _g_axes.append(_gene_panel(fig.add_subplot(gs[_eri]),GraphicRecord(sequence_length=_seqLen,features=_egfs),_eg['label'][:16],30,_is_last and not _has_cov))
            except Exception:
                _ax_eg=fig.add_subplot(gs[_eri]); _ax_eg.set_visible(False)
    if _HL_CODONS and _seq_sub:
        _all_hl=_fwd_hl+_rev_hl
        for _pax in _g_axes:
            _yl=_pax.get_ylim()
            for _x in _all_hl:
                _pax.axvspan(_x,_x+3,color='#3fb950',alpha=0.32,linewidth=0,zorder=3)
            _pax.set_ylim(_yl)
    if _FIG_CODON_COLS and _seq_sub:
        _u=_seq_sub.upper(); _fcd={}
        for _i in range(len(_u)-2):
            _cdn=_u[_i:_i+3]
            if _cdn in _FIG_CODON_COLS: _fcd.setdefault(_cdn,[]).append(_seq_sub_off+_i)
        for _cdn,_pxs in _fcd.items():
            _cc=_FIG_CODON_COLS[_cdn]
            for _pax in _g_axes:
                _yl=_pax.get_ylim()
                for _x in _pxs: _pax.axvline(_x,color=_cc,linewidth=0.8,alpha=0.7,zorder=3)
                _pax.set_ylim(_yl)
    if _CONTIG_SEPS:
        for _pax in _g_axes:
            for _sx in _CONTIG_SEPS:
                _pax.axvline(_sx,color='#58a6ff',linewidth=0.8,linestyle='--',alpha=0.6,zorder=4)
    buf=io.BytesIO()
    fig.savefig(buf,format='${fmt_}',dpi=${dpi},bbox_inches='tight',facecolor=_bg)
    plt.close(fig); buf.seek(0)
    _result='OK:'+base64.b64encode(buf.read()).decode()
except Exception:
    _result='ERR:'+traceback.format_exc()
`;
}

// ── Render ────────────────────────────────────────────────────────────────────
async function triggerRender(fmt) {
  const forPng = fmt==='png' || fmt===true;
  if (!pyReady||!lastState||pyRendering) return;
  pyRendering=true;
  document.getElementById('btn-render').disabled=true;
  setFigStatus(forPng?'rendering PNG\u2026':(fmt==='pdf'?'rendering PDF\u2026':'rendering\u2026'));
  const seq      = lastState.sequence || '';
  const seqStart = lastState.seq_start || 0;
  const relVS    = Math.max(0, localVS - seqStart);
  const relVE    = Math.min(seq.length, localVE - seqStart);
  const subseq   = relVS < relVE ? seq.slice(relVS, relVE) : '';
  const subOff   = Math.max(0, seqStart - localVS);
  pyodide.globals.set('_seq_sub', subseq);
  pyodide.globals.set('_seq_sub_off', subOff);
  try {
    await pyodide.runPythonAsync(buildPyCode(localVS,localVE,fmt));
    const result=pyodide.globals.get('_result');
    if (!result||typeof result!=='string') {
      setFig('<div class="fig-ph">Python returned null — check console.</div>');
      setFigStatus('error');
    } else if (result.startsWith('ERR:')) {
      setFig(`<pre class="fig-err">${esc(result.slice(4))}</pre>`);
      setFigStatus('Python error \u2014 see figure panel');
    } else if (result.startsWith('OK:')) {
      const b64=result.slice(3);
      if (forPng) {
        const blob=await (await fetch('data:image/png;base64,'+b64)).blob();
        dlBlob(blob,'genetui_figure.png'); setFigStatus('PNG downloaded');
      } else if (fmt==='pdf') {
        const blob=await (await fetch('data:application/pdf;base64,'+b64)).blob();
        dlBlob(blob,'genetui_figure.pdf'); setFigStatus('PDF downloaded');
      } else {
        lastSvgData=atob(b64);
        setFig(`<img src="data:image/svg+xml;base64,${b64}" style="max-width:100%">`);
        document.getElementById('btn-svg').disabled=false;
        document.getElementById('btn-png').disabled=false;
        document.getElementById('btn-pdf').disabled=false;
        document.getElementById('btn-py').disabled=false;
        setFigStatus('ready');
      }
    }
  } catch(err) {
    setFig(`<pre class="fig-err">${esc(String(err))}</pre>`);
    setFigStatus('JS error');
  }
  document.getElementById('btn-render').disabled=false;
  pyRendering=false;
}

// ── Button handlers ───────────────────────────────────────────────────────────
document.getElementById('btn-render').addEventListener('click',()=>triggerRender('svg'));
document.getElementById('btn-svg').addEventListener('click',()=>{
  if (lastSvgData) dlBlob(new Blob([lastSvgData],{type:'image/svg+xml'}),'genetui_figure.svg');
});
document.getElementById('btn-png').addEventListener('click',()=>triggerRender('png'));
document.getElementById('btn-pdf').addEventListener('click',()=>triggerRender('pdf'));
document.getElementById('btn-py').addEventListener('click',()=>{
  if (!lastState) return;
  const vs=localVS, ve=localVE;
  let code = buildPyCode(vs, ve, 'svg');
  if (!code) return;
  code = code
    .replace(/import io,base64,traceback,matplotlib\n/,'import traceback,matplotlib\n')
    .replace(/import io,base64,traceback\n/,'import traceback\n')
    .replace(/_result=None\n/,'')
    .replace(/try:\n/,'try:\n')
    .replace(/\s*_result='OK:'\+base64\.b64encode\(buf\.read\(\)\)\.decode\(\)\n/,
             "    plt.savefig('genetui_figure.png', dpi=600, bbox_inches='tight')\n    print('Saved genetui_figure.png')\n")
    .replace(/\s*_result='ERR:'\+traceback\.format_exc\(\)\n/,
             "    print(traceback.format_exc())\n")
    .replace(/buf\s*=\s*io\.BytesIO\(\)\n\s*[^\n]*\.savefig\(buf[^\n]*\)\n\s*buf\.seek\(0\)\n/g, '');
  const files = (lastState.input_files||[]).map(f=>`#   ${f}`).join('\n');
  const filesLine = files ? `# Source files:\n${files}\n` : '';
  const header = `# genetui figure — ${lastState.genome_name} ${vs}–${ve}\n# Generated by genetui (https://github.com/ZacharyArdern/genetui)\n${filesLine}# Run: pip install dna_features_viewer matplotlib && python genetui_figure.py\n\n`;
  dlBlob(new Blob([header+code],{type:'text/plain'}),'genetui_figure.py');
});
