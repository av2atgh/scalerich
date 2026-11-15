"""LeafVeinCNN per-leaf xlsx -> npz (bio/data/venation/<code>.npz). Run with .venv python -I."""
import glob, os, sys, numpy as np, openpyxl
src='/Users/avazquez/av2atg/scalerich/bio/data/raw/venation_leafveincnn/xlsx'
dst='/Users/avazquez/av2atg/scalerich/bio/data/venation'
def sheet(wb,name,cols):
    ws=wb[name]; rows=ws.iter_rows(values_only=True); h=list(next(rows)); idx=[h.index(c) for c in cols]
    out={c:[] for c in cols}
    for r in rows:
        if r[0] is None: continue
        for c,i in zip(cols,idx): out[c].append(r[i])
    return out
def num(x): return np.array([np.nan if v is None or isinstance(v,str) else v for v in x],float)
for f in sorted(glob.glob(src+'/*.xlsx')):
    code=os.path.basename(f).replace('_results.xlsx','')
    o=f'{dst}/{code}.npz'
    if os.path.exists(o): continue
    wb=openpyxl.load_workbook(f,read_only=True)
    E=sheet(wb,'Vein Edges',['EndNodes_1','EndNodes_2','Type','Length','Width','Or_ij','Or_ji','Ai','Aj','Tortuosity'])
    N=sheet(wb,'Vein Nodes',['node_ID','node_X_pix','node_Y_pix','node_Degree','node_D0','node_D1','node_D2','node_OD0','node_OD1','node_OD2'])
    H=sheet(wb,'HLD Nodes',['node_ID','width_threshold','node_Area','area_Asymmetry','Boundary','node_Degree'])
    HE=sheet(wb,'HLD Edges',['EndNodes_1','EndNodes_2'])
    P=sheet(wb,'Polygon Nodes',['polygon_ID','Area'])
    d={}
    for pre,S in (('e_',E),('n_',N),('h_',H),('he_',HE),('p_',P)):
        for c,v in S.items():
            d[pre+c]=np.array([1 if s=='EL' else 0 for s in v]) if c=='Type' else num(v)
    np.savez_compressed(o,**d); print(code,len(E['Length']),len(N['node_ID']),flush=True)
