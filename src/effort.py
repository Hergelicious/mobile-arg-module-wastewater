import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, json
from scipy.stats import spearmanr, kruskal
rng=np.random.default_rng(55)
x=pd.ExcelFile('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx')
s=pd.read_excel(x,'Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index); M=g.M_prior.values
def effort(sh):
    d=pd.read_excel(x,sh); sn=[c for c in d.columns if c.endswith('.sequence.number')]
    e=d[sn].max(); e.index=[c.replace('.sequence.number','') for c in e.index]; return e
eA=effort('Fig. S1a'); eC=effort('Fig. S1c')
# map rarefaction IDs to samples (Sample.ID.raw), verify 1:1
idmap=pd.Series(s.index,index=s.index)
E=pd.DataFrame({'effA':idmap.map(eA),'effC':idmap.map(eC)},index=s.index)
out={'n_matched_A':int(E.effA.notna().sum()),'n_matched_C':int(E.effC.notna().sum())}
G=s[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2)
pmap=s.WWTPID
LP=L.groupby(pmap).mean(); pm=s.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LP.index]
EP=np.log10(E).groupby(pmap).mean().loc[LP.index]
# is effort structured by country?
for c in ['effA','effC']:
    v=EP[c]; r2=1-((v-v.groupby(pm.country).transform('mean'))**2).sum()/((v-v.mean())**2).sum()
    out[f'{c}_country_R2']=float(r2)
    out[f'{c}_rho_mobile_share']=float(spearmanr(v, np.log10((G[g.index[g.M_prior]].sum(1)/G[core].sum(1)).groupby(pmap).mean().loc[LP.index]))[0])
def r2(Y,H):
    Yc=Y-Y.mean(0);P=H@np.linalg.pinv(H.T@H)@H.T;return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
dum=lambda l: pd.get_dummies(l).values.astype(float)
cc=pm.groupby('city').country.first()
def D_adj(Y):
    o=r2(Y,dum(pm.country)); d0=o[M].mean()-o[~M].mean(); nl=np.empty(999)
    for i in range(999):
        lab=pd.Series(rng.permutation(cc.values),index=cc.index); oo=r2(Y,dum(pm.city.map(lab))); nl[i]=oo[M].mean()-oo[~M].mean()
    return [float(d0),float((np.sum(nl>=d0)+1)/1000)]
def resid(Y,covs):
    X=np.c_[np.ones(len(covs)),covs]; B=np.linalg.lstsq(X,Y,rcond=None)[0]; return Y-X@B
Y=LP.values
out['D_unadjusted']=D_adj(Y)
for c in ['effA','effC']:
    out[f'D_adj_{c}']=D_adj(resid(Y,EP[[c]].values))
out['D_adj_both']=D_adj(resid(Y,EP[['effA','effC']].values))
json.dump(out,open('effort.json','w'),indent=1); print(json.dumps(out,indent=1))
