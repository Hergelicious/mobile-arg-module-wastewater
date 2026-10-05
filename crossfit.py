"""Cross-fitted module discovery and country-balanced discovery.
(1) CROSS-FITTING: in each of 200 iterations, 70% of cities are used to discover a module (within-country
    correlations, average linkage, k by silhouette, cluster most enriched for mobile-associated genes by
    hypergeometric test); membership is then FROZEN and evaluated in the held-out cities only.
    Held-out measures: mean within-country pairwise correlation of module genes vs non-module genes;
    D (mean country R2 of module minus non-module genes); reliability (ICC) of the module score.
(2) STABLE CORE: how often each gene is selected across iterations.
(3) COUNTRY-BALANCED: discovery repeated (a) with countries weighted equally when computing correlations and
    (b) excluding the USA and China entirely; agreement with the primary module reported as Jaccard overlap.
Nothing here modifies any previously reported result."""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.stats import hypergeom
from sklearn.metrics import silhouette_score
rng=np.random.default_rng(4242)
U='/mnt/user-data/uploads/'
s=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index); M=g.M_prior.values.astype(bool)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]; PRIM=set(mm.index[mm.in_module])
G=s[core].astype(float); pc=G.replace(0,np.nan).min()/2
L=np.log10(G+pc); LP=L.groupby(s.WWTPID).mean()
pm=s.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LP.index]
def discover(rows,weights_by_country=False):
    sub=LP.loc[rows]; ctry=pm.loc[rows,'country']
    if weights_by_country:
        Cs=[]
        for c,idx in sub.groupby(ctry.values).groups.items():
            d=sub.loc[idx]
            if len(d)>=4: Cs.append(np.corrcoef((d-d.mean()).values.T))
        C=np.nanmean(np.stack(Cs),0)
    else:
        d=sub-sub.groupby(ctry.values).transform('mean'); C=np.corrcoef(d.values.T)
    np.fill_diagonal(C,1.0); C=np.nan_to_num(C); D=np.clip(1-C,0,2)
    Z=linkage(squareform(D,checks=False),method='average')
    best=None
    for k in range(2,13):
        lab=fcluster(Z,k,criterion='maxclust')
        if len(set(lab))>1:
            sil=silhouette_score(D,lab,metric='precomputed')
            if best is None or sil>best[0]: best=(sil,lab)
    lab=best[1]; N=len(core); K=int(M.sum()); bestq=(1.1,None)
    for c in set(lab):
        idx=lab==c; n=int(idx.sum())
        if n>=5:
            p=hypergeom.sf(int((idx&M).sum())-1,N,K,n)
            if p<bestq[0]: bestq=(p,idx)
    return bestq[1], bestq[0]
def icc1(v,grp):
    d=pd.DataFrame({'v':v,'g':grp}); d=d[d.g.map(d.g.value_counts())>=2]
    k=d.g.value_counts(); n=len(k); N=len(d)
    msb=(k*(d.groupby('g').v.mean()-d.v.mean())**2).sum()/(n-1); msw=((d.v-d.groupby('g').v.transform('mean'))**2).sum()/(N-n)
    kbar=(N-(k**2).sum()/N)/(n-1); return float((msb-msw)/(msb+(kbar-1)*msw))
def r2(Y,H):
    Yc=Y-Y.mean(0); P=H@np.linalg.pinv(H.T@H)@H.T; return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
dum=lambda l: pd.get_dummies(pd.Series(l)).values.astype(float)
cities=pm.city.unique(); sel=np.zeros(len(core)); res=[]
for it in range(200):
    tr=set(rng.choice(cities,int(round(0.7*len(cities))),replace=False))
    rows_tr=[p for p in LP.index if pm.loc[p,'city'] in tr]; rows_te=[p for p in LP.index if pm.loc[p,'city'] not in tr]
    if len(rows_te)<20: continue
    mask,q=discover(rows_tr)
    if mask is None: continue
    sel+=mask
    te=LP.loc[rows_te]; cte=pm.loc[rows_te,'country']
    d=te-te.groupby(cte.values).transform('mean'); C=np.corrcoef(d.values.T); np.fill_diagonal(C,np.nan)
    iu=lambda m_: float(np.nanmean(C[np.ix_(m_,m_)]))
    o=r2(te.values,dum(cte.values))
    sc=L.loc[[i for i in s.index if s.WWTPID[i] in rows_te],[core[i] for i in np.where(mask)[0]]].mean(1)
    res.append({'k_module':int(mask.sum()),'q':float(q),'r_module':iu(mask),'r_other':iu(~mask),
                'D':float(o[mask].mean()-o[~mask].mean()),'ICC':icc1(sc,s.WWTPID[sc.index]),
                'jaccard_with_primary':len(set(np.array(core)[mask])&PRIM)/len(set(np.array(core)[mask])|PRIM)})
R=pd.DataFrame(res); freq=pd.Series(sel/len(R),index=core)
OUT={'n_iterations':int(len(R)),'held_out_summary':{k:[round(float(R[k].median()),3),round(float(R[k].quantile(.05)),3),round(float(R[k].quantile(.95)),3)] for k in ['k_module','r_module','r_other','D','ICC','jaccard_with_primary']},
     'D_positive_fraction':float((R.D>0).mean()),
     'stable_core_90':sorted(freq[freq>=0.9].index),'stable_core_75':sorted(freq[freq>=0.75].index),
     'selection_frequency':{k:round(float(v),3) for k,v in freq.sort_values(ascending=False).items()}}
# country-balanced variants
mask_w,_=discover(list(LP.index),weights_by_country=True)
rows_ex=[p for p in LP.index if pm.loc[p,'country'] not in ('USA','China')]
mask_ex,_=discover(rows_ex)
J=lambda a,b: len(a&b)/len(a|b)
OUT['country_balanced']={'equal_weight_module_size':int(mask_w.sum()),'equal_weight_jaccard':round(J(set(np.array(core)[mask_w]),PRIM),3),
  'equal_weight_mobile_captured':int((mask_w&M).sum()),
  'excluding_USA_China_module_size':int(mask_ex.sum()),'excluding_USA_China_jaccard':round(J(set(np.array(core)[mask_ex]),PRIM),3),
  'excluding_USA_China_mobile_captured':int((mask_ex&M).sum()),'n_plants_excluding':len(rows_ex)}
json.dump(OUT,open('crossfit.json','w'),indent=1)
print(json.dumps({k:v for k,v in OUT.items() if k!='selection_frequency'},indent=1)[:1500])
