import json, hashlib, numpy as np, pandas as pd
from scipy.stats import hypergeom
txt=open('frozen_spec2.json').read()
assert hashlib.sha256(txt.encode()).hexdigest()=='a2f3b77ed1f2e46b823cc82714dd05b65dd0007e91e660639c154b9b22486be0','spec changed!'
S=json.loads(txt)['D_external_classification']
rng=np.random.default_rng(90210)
s=pd.read_excel('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]
RF=np.array([x in set(S['resfinder_acquired_genes']) for x in core]); CUR=g.M_prior.values.astype(bool); MOD=mm.in_module.values.astype(bool)
assert RF.sum()==S['n_acquired']
G=s[core].astype(float); LP=np.log10(G+G.replace(0,np.nan).min()/2).groupby(s.WWTPID).mean()
pm=s.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LP.index]
def r2(Y,H):
    Yc=Y-Y.mean(0); P=H@np.linalg.pinv(H.T@H)@H.T; return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
dum=lambda l: pd.get_dummies(l).values.astype(float)
Y=LP.values; obs=r2(Y,dum(pm.country))
K=json.load(open('canonical.json')); assert abs((obs[CUR].mean()-obs[~CUR].mean())-K['D'][0])<1e-12  # same matrix as primary
cc=pm.groupby('city').country.first(); NP=9999; N=np.empty((NP,len(core)))
for i in range(NP):
    lab=pd.Series(rng.permutation(cc.values),index=cc.index); N[i]=r2(Y,dum(pm.city.map(lab)))
def test(a,b):
    o=obs[a].mean()-obs[b].mean(); n=N[:,a].mean(1)-N[:,b].mean(1); return [float(o),float((np.sum(n>=o)+1)/(NP+1))]
OUT={'spec_sha256':hashlib.sha256(txt.encode()).hexdigest(),'n_acquired':int(RF.sum()),'n_other':int((~RF).sum())}
OUT['primary_D_RF']=test(RF,~RF)
OUT['mean_R2']={'resfinder':float(obs[RF].mean()),'other':float(obs[~RF].mean())}
# secondary 1: agreement with curated list
a=int((RF&CUR).sum()); b=int((RF&~CUR).sum()); c=int((~RF&CUR).sum()); d=int((~RF&~CUR).sum()); n=a+b+c+d
po=(a+d)/n; pe=((a+b)*(a+c)+(c+d)*(b+d))/n**2
OUT['agreement']={'both':a,'resfinder_only':b,'curated_only':c,'neither':d,'kappa':float((po-pe)/(1-pe))}
OUT['curated_only_genes']=[core[i] for i in np.where(~RF&CUR)[0]]; OUT['resfinder_only_genes']=[core[i] for i in np.where(RF&~CUR)[0]]
# secondary 2: beyond our labels (among the 69 non-curated genes)
OUT['D_RF_noncurated']=test(RF&~CUR,~RF&~CUR); OUT['n_noncurated_split']=[int((RF&~CUR).sum()),int((~RF&~CUR).sum())]
# secondary 3: module enrichment
Nn,Kk,nn,kk=len(core),int(RF.sum()),int(MOD.sum()),int((RF&MOD).sum())
OUT['module_enrichment']={'module_size':nn,'resfinder_in_module':kk,'P':float(hypergeom.sf(kk-1,Nn,Kk,nn))}
OUT['success']=bool(OUT['primary_D_RF'][1]<0.05 and OUT['primary_D_RF'][0]>0)
json.dump(OUT,open('frozen_results2.json','w'),indent=1); print(json.dumps(OUT,indent=1))
