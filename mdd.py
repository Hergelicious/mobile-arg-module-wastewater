"""Minimum detectable difference for the module score, from the measured variance components.
Discovery cohort: within-plant (replicate) SD and between-plant SD of the module score.
Reports the smallest difference detectable between two plants, and the smallest change within one plant between
two occasions, at 80% power and alpha=0.05, for 1-5 samples per plant; in log10 units and as a fold change."""
import warnings; warnings.filterwarnings("ignore")
import json, numpy as np, pandas as pd
U='/mnt/user-data/uploads/'
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]
MOD=[x for x in core if bool(mm.loc[x,'in_module'])]
z=pd.read_excel(U+'41467_2025_59019_MOESM8_ESM.xlsx','Fig. S2b and Fig. S4',index_col=0)
G=z[core].astype(float); L=np.log10(G+G.replace(0,np.nan).min()/2)
score=L[MOD].mean(1)
d=pd.DataFrame({'v':score,'p':z.WWTPID})
rep=d[d.p.map(d.p.value_counts())>=2]
k=rep.p.value_counts(); n=len(k); N=len(rep)
msw=((rep.v-rep.groupby('p').v.transform('mean'))**2).sum()/(N-n)
sd_w=float(np.sqrt(msw)); sd_b=float(d.groupby('p').v.mean().std())
OUT={'within_plant_SD_log10':round(sd_w,3),'between_plant_SD_log10':round(sd_b,3),'n_plants_with_replicates':int(n),
     'module_score_range_log10':round(float(d.groupby('p').v.mean().max()-d.groupby('p').v.mean().min()),2),'MDD':{}}
for m_ in [1,2,3,5]:
    se=sd_w*np.sqrt(2.0/m_)                     # difference of two plant means, m samples each
    mdd=2.8*se                                   # 80% power, two-sided alpha 0.05
    OUT['MDD'][m_]={'log10':round(float(mdd),3),'fold':round(float(10**mdd),2),
                    'as_%_of_observed_range':round(100*float(mdd)/OUT['module_score_range_log10'],1)}
json.dump(OUT,open('mdd.json','w'),indent=1); print(json.dumps(OUT,indent=1))
