import warnings; warnings.filterwarnings("ignore")
import json, hashlib, itertools, numpy as np, pandas as pd
from scipy.stats import spearmanr, rankdata
txt=open('frozen_spec.json').read()
assert hashlib.sha256(txt.encode()).hexdigest()=='57666824bf6f6c7f14a87bc989f290a921fa0e60cfd5a70c57df0eaaf7754b7d', 'spec changed!'
S=json.loads(txt); OUT={'spec_sha256':hashlib.sha256(txt.encode()).hexdigest()}
rng=np.random.default_rng(20260920)
# ---------------- A: network co-membership ----------------
A=S['A_network_test']
nodes=pd.read_excel('/mnt/user-data/uploads/41467_2025_66070_MOESM6_ESM.xlsx')
acq=nodes[(nodes.source=='amr')&(nodes.type=='acquired')].set_index('taxa_id')
mod_nodes=sorted({A['mapping'][gname] for gname in A['module_genes_mapped']})
k=len(mod_nodes); OUT['A_n_module_nodes']=k
if k<6: OUT['A_result']='uninformative (fewer than 6 module nodes)'
else:
    comm=acq.comm_id
    def share(ns):
        pairs=list(itertools.combinations(ns,2)); return np.mean([comm[a]==comm[b] for a,b in pairs])
    links=pd.read_excel('/mnt/user-data/uploads/41467_2025_66070_MOESM7_ESM.xlsx',usecols=['from','to','link_type'])
    aa=links[links.link_type=='acquired-acquired']; E={frozenset(p) for p in zip(aa['from'],aa['to'])}
    def density(ns):
        pairs=list(itertools.combinations(ns,2)); return np.mean([frozenset(p) in E for p in pairs])
    obs=share(mod_nodes); obsd=density(mod_nodes)
    pool=list(acq.index); nullS=np.empty(10000); nullD=np.empty(10000)
    for i in range(10000):
        pick=list(rng.choice(pool,k,replace=False)); nullS[i]=share(pick); nullD[i]=density(pick)
    OUT['A_primary']={'observed_share':float(obs),'null_mean':float(nullS.mean()),'null_q95':float(np.quantile(nullS,.95)),'P':float((np.sum(nullS>=obs)+1)/10001)}
    OUT['A_secondary']={'observed_density':float(obsd),'null_mean':float(nullD.mean()),'P':float((np.sum(nullD>=obsd)+1)/10001)}
    OUT['A_n_communities_among_module_nodes']=int(comm[mod_nodes].nunique()); OUT['A_n_communities_among_acquired']=int(comm.nunique())
    OUT['A_success']=bool(OUT['A_primary']['P']<0.05)
# ---------------- B: intI1 ----------------
x=pd.ExcelFile('/mnt/user-data/uploads/41467_2025_59019_MOESM5_ESM.xlsx')
mge=pd.read_excel(x,header=2); mge=mge.rename(columns={mge.columns[0]:'orf'})
rows=mge[mge['MGE_Gene_name'].astype(str).str.lower()=='inti1']
OUT['B_n_intI1_rows']=int(len(rows))
s=pd.read_excel('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index); M=g.M_prior.values
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]; MOD=mm.in_module.values
if len(rows):
    samp=[c for c in rows.columns if c in s.index]
    OUT['B_samples_matched']=len(samp)
    intI=rows[samp].astype(float).sum(0).reindex(s.index)
    G=s[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2)
    li=np.log10(intI+intI[intI>0].min()/2)
    LP=L.groupby(s.WWTPID).mean(); IP=li.groupby(s.WWTPID).mean().loc[LP.index]
    ctry=s.groupby('WWTPID').CountryRegion_full.first().loc[LP.index]
    Ld=LP-LP.groupby(ctry.values).transform('mean'); Id=IP-IP.groupby(ctry.values).transform('mean')
    RL=np.apply_along_axis(rankdata,0,Ld.values)
    def rhos(iv):
        ri=rankdata(iv); ri=(ri-ri.mean())/ri.std(); Z=(RL-RL.mean(0))/RL.std(0); return (Z*ri[:,None]).mean(0)
    r=rhos(Id.values); stat=lambda rr,mask: rr[mask].mean()-rr[~mask].mean()
    o1=stat(r,M); o2=stat(r,MOD)
    grp=ctry.values; idx_by=[np.where(grp==c)[0] for c in np.unique(grp)]
    n1=np.empty(9999); n2=np.empty(9999); iv=Id.values
    for i in range(9999):
        perm=iv.copy()
        for ix in idx_by: perm[ix]=iv[rng.permutation(ix)]
        rr=rhos(perm); n1[i]=stat(rr,M); n2[i]=stat(rr,MOD)
    OUT['B_primary']={'mean_rho_mobile':float(r[M].mean()),'mean_rho_other':float(r[~M].mean()),'difference':float(o1),'P':float((np.sum(n1>=o1)+1)/10000)}
    OUT['B_secondary']={'mean_rho_module':float(r[MOD].mean()),'mean_rho_outside':float(r[~MOD].mean()),'difference':float(o2),'P':float((np.sum(n2>=o2)+1)/10000)}
    OUT['B_success']=bool(OUT['B_primary']['P']<0.05 and o1>0)
# ---------------- C: depth adjustment ----------------
d1=pd.read_excel('/mnt/user-data/uploads/41467_2025_59019_MOESM3_ESM.xlsx',header=2).set_index('Sample')
OUT['C_samples_matched']=int(d1.index.isin(s.index).sum())
hq=np.log10(d1['HQ bases'].reindex(s.index).astype(float))
G=s[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2); LP=L.groupby(s.WWTPID).mean()
HP=hq.groupby(s.WWTPID).mean().loc[LP.index]; pm=s.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LP.index]
Y=LP.values; X=np.c_[np.ones(len(HP)),HP.values]; Yr=Y-X@np.linalg.lstsq(X,Y,rcond=None)[0]
def r2(Yv,H):
    Yc=Yv-Yv.mean(0); P=H@np.linalg.pinv(H.T@H)@H.T; return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
dum=lambda l: pd.get_dummies(l).values.astype(float); cc=pm.groupby('city').country.first()
o=r2(Yr,dum(pm.country)); d0=o[M].mean()-o[~M].mean(); nl=np.empty(999)
for i in range(999):
    lab=pd.Series(rng.permutation(cc.values),index=cc.index); oo=r2(Yr,dum(pm.city.map(lab))); nl[i]=oo[M].mean()-oo[~M].mean()
hqc=HP.groupby(pm.country.values).mean()
OUT['C_depth']={'HQ_Gb_median':float(10**hq.median()/1e9),'HQ_Gb_range':[float(10**hq.min()/1e9),float(10**hq.max()/1e9)],
   'depth_country_R2':float(1-((HP-HP.groupby(pm.country.values).transform('mean'))**2).sum()/((HP-HP.mean())**2).sum()),
   'D_depth_adjusted':float(d0),'P':float((np.sum(nl>=d0)+1)/1000)}
json.dump(OUT,open('frozen_results.json','w'),indent=1); print(json.dumps(OUT,indent=1))
