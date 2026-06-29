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
function buildPyCode(vs,ve,forPng) {
  if (!lastState) return null;
  const width     = parseFloat(document.getElementById('fig-width').value);
  const fs        = parseInt(document.getElementById('fig-fs').value);
  const showLbls  = document.getElementById('fig-show-labels').checked;
  const hideLong  = document.getElementById('fig-hide-long-labels').checked;
  const ruler     = document.getElementById('fig-ruler').checked;
  const white     = document.getElementById('fig-white').checked;
  const sixframe  = document.getElementById('opt-sixframe').checked;
  const stopCod   = document.getElementById('opt-stopcodons').checked;
  const seqLen    = Math.max(1,ve-vs);
  const feats     = lastState.features.filter(f=>f.start<=ve&&f.end>=vs);
  const bg        = white?"'white'":"'#0d1117'";
  const fg        = white?"'#111'":"'#c9d1d9'";
  const dpi       = forPng?600:96;
  const fmt_      = forPng?'png':'svg';

  const covStyle   = document.getElementById('cov-style-fig').value;
  const covBinSz   = lastState.coverage_bin_size  || 1000;
  const covBinOff  = lastState.coverage_bin_start || 0;
  const covPlus    = lastState.coverage_plus  || [];
  const covMinus   = lastState.coverage_minus || [];
  let pyCovPlus = '[]', pyCovMinus = '[]', pyCovBinSz = covBinSz, pyCovBinOff = 0;
  if (covStyle !== 'none' && covPlus.length) {
    const fb = Math.max(0, Math.floor((vs - covBinOff) / covBinSz));
    const lb = Math.min(covPlus.length - 1, Math.ceil((ve - covBinOff) / covBinSz));
    if (fb <= lb) {
      pyCovPlus   = JSON.stringify(covPlus.slice(fb, lb+1));
      pyCovMinus  = JSON.stringify(covMinus.slice(fb, lb+1));
      pyCovBinOff = covBinOff + fb * covBinSz;
    }
  }
  const hasCov = covStyle !== 'none' && pyCovPlus !== '[]';

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
    _fwd_stops=[[],[],[]]
    _rev_stops=[[],[],[]]
    if _stopCod and _seq_sub:
        _u=_seq_sub.upper()
        _rc_u=_rc(_seq_sub)
        for _i in range(len(_u)-2):
            if _u[_i:_i+3] in _STOP:
                _gfr=(_geno_off+_i)%3
                _fwd_stops[_gfr].append(_seq_sub_off+_i)
        for _i in range(len(_rc_u)-2):
            if _rc_u[_i:_i+3] in _STOP:
                _gend=_geno_off+len(_seq_sub)-_i-1
                _gfr=_gend%3
                _rev_stops[_gfr].append(_seq_sub_off+len(_seq_sub)-_i-2)
    _ROWS=[('+1',1,0),('+2',1,1),('+3',1,2),('-1',-1,0),('-2',-1,1),('-3',-1,2)]
    _cov_style='${covStyle}'; _cov_plus=${pyCovPlus}; _cov_minus=${pyCovMinus}
    _cov_bin_sz=${pyCovBinSz}; _cov_bin_off=${pyCovBinOff}
    _has_cov=_cov_style!='none' and len(_cov_plus)>0
    # row order: cov+? / nc+ / +1/+2/+3 / sep / -1/-2/-3 / nc- / cov-?
    _hr=([1.4] if _has_cov else [])+[0.45,0.75,0.75,0.75,0.25,0.75,0.75,0.75,0.45]+([1.4] if _has_cov else [])
    _off=1 if _has_cov else 0
    _fig_h=6*0.29+2*0.22+0.75+(1.6 if _has_cov else 0)
    fig=plt.figure(figsize=(_width,_fig_h))
    fig.patch.set_facecolor(_bg)
    gs=gridspec.GridSpec(len(_hr),1,figure=fig,hspace=0.0,
        height_ratios=_hr,top=0.96,bottom=0.09,left=0.11,right=0.98)
    import matplotlib.ticker as _tck
    def _gfmt(x,p):
        v=int(x)+${vs}
        if v>=1000000: return f'{v//1000000}M'
        if v>=1000: return f'{v//1000}k'
        return str(v)
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
            _gfs.append(GraphicFeature(start=s,end=e,strand=_si,label=lb,color=col))
        rec=GraphicRecord(sequence_length=_seqLen,features=_gfs)
        rec.plot(ax=ax,with_ruler=False,draw_line=True)
        ax.set_ylim(-0.55,0.55)
        ax.set_facecolor(_bg)
        ax.set_ylabel(_lbl,fontsize=7,rotation=0,labelpad=22,va='center',color=_fg)
        for sp in ax.spines.values(): sp.set_color(_fg)
        ax.tick_params(colors=_fg,labelsize=${fs-2})
        ax.set_xticklabels([])
        if _stopCod and _seq_sub:
            _sc_list=_fwd_stops[_fr] if _si==1 else _rev_stops[_fr]
            for _x in _sc_list:
                ax.axvline(_x,color='black',alpha=0.7,linewidth=0.8)
        _axes.append(ax)
    _nc_p=[GraphicFeature(start=s,end=e,strand=st,label=lb,color=col) for (s,e,st,lb,col,nc,fr) in _feats_raw if nc and st==1]
    _nc_m=[GraphicFeature(start=s,end=e,strand=st,label=lb,color=col) for (s,e,st,lb,col,nc,fr) in _feats_raw if nc and st==-1]
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
        _xs=[(_cov_bin_off+i*_cov_bin_sz-${vs}) for i in range(len(_cov_plus))]
        _lp=_np.log10(_np.array(_cov_plus,dtype=float)+1)
        _lm=_np.log10(_np.array(_cov_minus,dtype=float)+1)
        _mxp=max(float(_lp.max()) if len(_lp) else 1,0.01)
        _mxm=max(float(_lm.max()) if len(_lm) else 1,0.01)
        def _cov_panel(ax,xs,raw,col,label,invert=False,show_xaxis=False):
            ax.set_facecolor(_bg)
            for sp in ax.spines.values(): sp.set_color(_fg)
            ax.tick_params(colors=_fg,labelsize=${fs-2})
            ax.set_ylabel(label,fontsize=7,rotation=0,labelpad=38,va='center',color=_fg)
            ax.set_xlim(0,${seqLen})
            if show_xaxis:
                ax.xaxis.set_major_formatter(_bpfmt)
                ax.tick_params(axis='x',colors=_fg,labelsize=${fs-2})
            else:
                ax.set_xticklabels([])
            if _cov_style=='reads':
                from matplotlib.patches import Rectangle
                from matplotlib.collections import PatchCollection
                _rh,_rg,_mr=0.8,0.3,20
                def _lcg(s): return (1664525*s+1013904223)&0xFFFFFFFF
                _rects=[]
                for _i,(_x,_c) in enumerate(zip(xs,raw)):
                    if _c==0: continue
                    _sd=_i*(67890 if invert else 12345)
                    for _r in range(min(_mr,int(_c))):
                        _sd=_lcg(_sd)
                        _rw=max(1.0,_cov_bin_sz*(0.55+((_sd&0xff)/256-0.5)*0.3))
                        _rx=_x+max(0,((_sd>>8)&0xff)/256*(_cov_bin_sz-_rw))
                        _rects.append(Rectangle((_rx,_r*(_rh+_rg)),_rw,_rh))
                if _rects: ax.add_collection(PatchCollection(_rects,facecolor=col,alpha=0.65,linewidth=0))
                ax.set_ylim(0,_mr*(_rh+_rg)); ax.autoscale_view()
                _rd_ticks=[v for v in [5,10,20] if v<=_mr]
                ax.set_yticks([v*(_rh+_rg) for v in _rd_ticks])
                ax.set_yticklabels([str(v) for v in _rd_ticks],fontsize=6)
            else:
                lv=_np.log10(_np.array(raw,dtype=float)+1)
                mx=max(float(lv.max()) if len(lv) else 1,0.01)
                if _cov_style=='histogram':
                    ax.bar(xs,list(lv),width=_cov_bin_sz*0.9,align='edge',color=col,alpha=0.75)
                else:
                    _sigma=0.8; _r=int(_np.ceil(3*_sigma))
                    _k=_np.exp(-0.5*(_np.arange(-_r,_r+1)/_sigma)**2); _k/=_k.sum()
                    _sm=_np.convolve(lv,_k,mode='same')
                    ax.fill_between(xs,_sm,alpha=0.2,color=col); ax.plot(xs,_sm,color=col,linewidth=1.2)
                ax.set_ylim(0,mx*1.1)
                _ytv=[v for v in [1,10,100,1000,10000] if _np.log10(v+1)<=mx*1.05]
                ax.set_yticks([_np.log10(v+1) for v in _ytv])
                ax.set_yticklabels([str(v) for v in _ytv],fontsize=6)
            if invert: ax.invert_yaxis()
        _cov_panel(fig.add_subplot(gs[0]),_xs,list(_np.array(_cov_plus)),'#58a6ff','cov+')
        _cov_panel(fig.add_subplot(gs[-1]),_xs,list(_np.array(_cov_minus)),'#f78166','cov-',invert=True,show_xaxis=True)
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
    _f_c   =[GraphicFeature(start=s,end=e,strand=si,label=lb,color=col) for (s,e,si,lb,col,nc,isp) in _feats_raw if not nc]
    _f_nc_p=[GraphicFeature(start=s,end=e,strand=si,label=lb,color=col) for (s,e,si,lb,col,nc,isp) in _feats_raw if nc and isp]
    _f_nc_m=[GraphicFeature(start=s,end=e,strand=si,label=lb,color=col) for (s,e,si,lb,col,nc,isp) in _feats_raw if nc and not isp]
    _has_c=bool(_f_c); _has_nc_p=bool(_f_nc_p); _has_nc_m=bool(_f_nc_m)
    def _gfmt(x,p):
        v=int(x)+${vs}
        if v>=1000000: return f'{v//1000000}M'
        if v>=1000: return f'{v//1000}k'
        return str(v)
    _bpfmt=_tck.FuncFormatter(_gfmt)
    _cov_style='${covStyle}'; _cov_plus=${pyCovPlus}; _cov_minus=${pyCovMinus}
    _cov_bin_sz=${pyCovBinSz}; _cov_bin_off=${pyCovBinOff}
    _has_cov=_cov_style!='none' and len(_cov_plus)>0
    # Assign explicit gridspec row indices: cov+ / nc+ / coding / nc- / cov-
    _nrows=0; _hr=[]
    _ri_cov_p=_ri_nc_p=_ri_c=_ri_nc_m=_ri_cov_m=-1
    if _has_cov:  _ri_cov_p=_nrows; _nrows+=1; _hr.append(1.4)
    if _has_nc_p: _ri_nc_p =_nrows; _nrows+=1; _hr.append(0.65)
    if _has_c:    _ri_c    =_nrows; _nrows+=1; _hr.append(2.5)
    if _has_nc_m: _ri_nc_m =_nrows; _nrows+=1; _hr.append(0.65)
    if _has_cov:  _ri_cov_m=_nrows; _nrows+=1; _hr.append(1.4)
    _fig_h=sum(_hr)+0.55
    fig=plt.figure(figsize=(_width,_fig_h)); fig.patch.set_facecolor(_bg)
    gs=gridspec.GridSpec(_nrows,1,figure=fig,height_ratios=_hr,hspace=0.04,top=0.97,bottom=0.09,left=0.11,right=0.98)
    import numpy as _np
    _xs=[(_cov_bin_off+i*_cov_bin_sz-${vs}) for i in range(len(_cov_plus))] if _has_cov else []
    def _cov_panel(ax,raw,col,label,invert=False,show_xaxis=False):
        ax.set_facecolor(_bg)
        for sp in ax.spines.values(): sp.set_color(_fg)
        ax.tick_params(colors=_fg,labelsize=${fs-2})
        ax.set_ylabel(label,fontsize=7,rotation=0,labelpad=38,va='center',color=_fg)
        ax.set_xlim(0,_seqLen)
        if _cov_style=='reads':
            from matplotlib.patches import Rectangle; from matplotlib.collections import PatchCollection
            _rh,_rg,_mr=0.8,0.3,20; _lcg=lambda s:(1664525*s+1013904223)&0xFFFFFFFF; _rects=[]
            for _i,(_x,_c) in enumerate(zip(_xs,raw)):
                if _c==0: continue
                _sd=_i*(67890 if invert else 12345)
                for _r2 in range(min(_mr,int(_c))):
                    _sd=_lcg(_sd); _rw=max(1.0,_cov_bin_sz*(0.55+((_sd&0xff)/256-0.5)*0.3))
                    _rx=_x+max(0,((_sd>>8)&0xff)/256*(_cov_bin_sz-_rw)); _rects.append(Rectangle((_rx,_r2*(_rh+_rg)),_rw,_rh))
            if _rects: ax.add_collection(PatchCollection(_rects,facecolor=col,alpha=0.65,linewidth=0))
            ax.set_ylim(0,_mr*(_rh+_rg)); ax.autoscale_view()
        else:
            lv=_np.log10(_np.array(raw,dtype=float)+1); mx=max(float(lv.max()) if len(lv) else 1,0.01)
            if _cov_style=='histogram': ax.bar(_xs,list(lv),width=_cov_bin_sz*0.9,align='edge',color=col,alpha=0.75)
            else:
                _sigma=0.8; _r=int(_np.ceil(3*_sigma)); _k=_np.exp(-0.5*(_np.arange(-_r,_r+1)/_sigma)**2); _k/=_k.sum()
                _sm=_np.convolve(lv,_k,mode='same'); ax.fill_between(_xs,_sm,alpha=0.2,color=col); ax.plot(_xs,_sm,color=col,linewidth=1.2)
            ax.set_ylim(0,mx*1.1)
            _ytv=[v for v in [1,10,100,1000,10000] if _np.log10(v+1)<=mx*1.05]
            ax.set_yticks([_np.log10(v+1) for v in _ytv]); ax.set_yticklabels([str(v) for v in _ytv],fontsize=6)
        if invert: ax.invert_yaxis()
        if show_xaxis:
            ax.xaxis.set_major_locator(_tck.AutoLocator()); ax.xaxis.set_major_formatter(_bpfmt)
            ax.tick_params(axis='x',colors=_fg,labelsize=${fs-2})
        else:
            ax.set_xticklabels([])
    def _gene_panel(ax,rec,label,lpad,show_xaxis):
        # with_ruler=True keeps x-axis enabled; with_ruler=False calls ax.axis('off')
        rec.plot(ax=ax,with_ruler=show_xaxis,draw_line=True)
        ax.set_facecolor(_bg)
        for sp in ax.spines.values(): sp.set_color(_fg)
        ax.set_ylabel(label,fontsize=7,rotation=0,labelpad=lpad,va='center',color=_fg)
        for t in ax.texts: t.set_fontsize(${fs})
        if show_xaxis:
            # override dna_features_viewer's relative-coord formatter with absolute-coord one
            ax.xaxis.set_major_formatter(_bpfmt)
            ax.tick_params(axis='x',colors=_fg,labelsize=${fs-2})
    _nc_p_show_x =_has_nc_p and not _has_c and not _has_nc_m and not _has_cov
    _c_show_x    =_has_c    and not _has_nc_m and not _has_cov
    _nc_m_show_x =_has_nc_m and not _has_cov
    if _ri_cov_p>=0: _cov_panel(fig.add_subplot(gs[_ri_cov_p]),_cov_plus, '#58a6ff','cov+')
    if _ri_nc_p >=0: _gene_panel(fig.add_subplot(gs[_ri_nc_p]), GraphicRecord(sequence_length=_seqLen,features=_f_nc_p),'nc+',22,_nc_p_show_x)
    if _ri_c    >=0: _gene_panel(fig.add_subplot(gs[_ri_c]),    GraphicRecord(sequence_length=_seqLen,features=_f_c),    'genes',30,_c_show_x)
    if _ri_nc_m >=0: _gene_panel(fig.add_subplot(gs[_ri_nc_m]), GraphicRecord(sequence_length=_seqLen,features=_f_nc_m),'nc-',22,_nc_m_show_x)
    if _ri_cov_m>=0: _cov_panel(fig.add_subplot(gs[_ri_cov_m]),_cov_minus,'#f78166','cov-',invert=True,show_xaxis=True)
    buf=io.BytesIO()
    fig.savefig(buf,format='${fmt_}',dpi=${dpi},bbox_inches='tight',facecolor=_bg)
    plt.close(fig); buf.seek(0)
    _result='OK:'+base64.b64encode(buf.read()).decode()
except Exception:
    _result='ERR:'+traceback.format_exc()
`;
}

// ── Render ────────────────────────────────────────────────────────────────────
async function triggerRender(forPng) {
  if (!pyReady||!lastState||pyRendering) return;
  pyRendering=true;
  document.getElementById('btn-render').disabled=true;
  setFigStatus(forPng?'rendering PNG\u2026':'rendering\u2026');
  const seq      = lastState.sequence || '';
  const seqStart = lastState.seq_start || 0;
  const relVS    = Math.max(0, localVS - seqStart);
  const relVE    = Math.min(seq.length, localVE - seqStart);
  const subseq   = relVS < relVE ? seq.slice(relVS, relVE) : '';
  const subOff   = Math.max(0, seqStart - localVS);
  pyodide.globals.set('_seq_sub', subseq);
  pyodide.globals.set('_seq_sub_off', subOff);
  try {
    await pyodide.runPythonAsync(buildPyCode(localVS,localVE,forPng));
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
      } else {
        lastSvgData=atob(b64);
        setFig(`<img src="data:image/svg+xml;base64,${b64}" style="max-width:100%">`);
        document.getElementById('btn-svg').disabled=false;
        document.getElementById('btn-png').disabled=false;
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
document.getElementById('btn-render').addEventListener('click',()=>triggerRender(false));
document.getElementById('btn-svg').addEventListener('click',()=>{
  if (lastSvgData) dlBlob(new Blob([lastSvgData],{type:'image/svg+xml'}),'genetui_figure.svg');
});
document.getElementById('btn-png').addEventListener('click',()=>triggerRender(true));
document.getElementById('btn-py').addEventListener('click',()=>{
  if (!lastState) return;
  const vs=localVS, ve=localVE;
  let code = buildPyCode(vs, ve, false);
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
