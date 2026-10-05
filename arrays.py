import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr
rng=np.random.default_rng(20260919)
x=pd.ExcelFile('/mnt/user-data/uploads/41467_2025_59019_MOESM8_ESM.xlsx')
s=pd.read_excel(x,'Fig. S2b and Fig. S4',index_col=0)
g=pd.read_csv('gene_table2.csv',index_col=0); core=list(g.index); M=g.M_prior.values
G=s[core].astype(float); LP=np.log10(G+G.replace(0,np.nan).min()/2).groupby(s.WWTPID).mean()
pm=s.groupby('WWTPID').agg(country=('CountryRegion_full','first'),city=('City','first')).loc[LP.index]
Y=LP.values; cc=pm.groupby('city').country.first()
def r2(Yv,H):
    Yc=Yv-Yv.mean(0);P=H@np.linalg.pinv(H.T@H)@H.T;return 1-((Yc-P@Yc)**2).sum(0)/(Yc**2).sum(0)
dum=lambda l: pd.get_dummies(l).values.astype(float)
obs=r2(Y,dum(pm.country)); N=np.empty((9999,len(core)))
for i in range(9999):
    lab=pd.Series(rng.permutation(cc.values),index=cc.index); N[i]=r2(Y,dum(pm.city.map(lab)))
np.save('Dnull.npy',N[:,M].mean(1)-N[:,~M].mean(1))
exc=obs-N.mean(0); np.save('excess.npy',exc)
cpl=g.M_coupling.values
np.save('rhonull.npy',np.array([spearmanr(cpl,N[i]-N.mean(0))[0] for i in range(2000)]))
g['excess_R2']=exc; g['R2_country']=obs
pg=((N>=obs).sum(0)+1)/10000; g['perm_p']=pg
g.to_csv('gene_table_final.csv'); print('saved', float((obs[M].mean()-obs[~M].mean())))
