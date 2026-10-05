"""(D) Design guidance for future wastewater resistome surveys (calibrated simulation; see header of benchmark for scope).
Model: y = country + city(country) + plant + season. Variance shares from this study's decomposition;
season share from the external cohort (27% of module-score variance).
Designs: 'aliased' (each country sampled in one window, as here) vs 'crossed' (windows spread across countries).
Reports power to detect the mobile-vs-other contrast and the country R2 estimate under each design."""
import warnings; warnings.filterwarnings("ignore")
import json, sys, numpy as np
rng=np.random.default_rng(31337)
K=json.load(open('canonical.json'))
sh={'mob':np.array([K['SS_country'][0],K['SS_city'][0],1-K['SS_country'][0]-K['SS_city'][0]]),
    'oth':np.array([K['SS_country'][1],K['SS_city'][1],1-K['SS_country'][1]-K['SS_city'][1]])}
SEASON=json.load(open('frozen_results3_alr.json'))['E_temporal']['variance_fractions']['year']
NM,NO=45,69; M=np.array([True]*NM+[False]*NO)
def r2_fast(Y,lab,ng):
    cnt=np.bincount(lab,minlength=ng).astype(float)
    s=np.zeros((ng,Y.shape[1]))
    np.add.at(s,lab,Y)
    mu=s/cnt[:,None]
    resid=Y-mu[lab]
    Yc=Y-Y.mean(0)
    return 1-(resid**2).sum(0)/(Yc**2).sum(0)
def simulate(k,cities,plants=2,design='aliased',nsim=60,nperm=49):
    ctry=np.repeat(np.arange(k),cities*plants); city=np.repeat(np.arange(k*cities),plants); n=len(ctry)
    season_key=ctry if design=='aliased' else np.tile(np.arange(4),int(np.ceil(n/4)))[:n]
    ns=season_key.max()+1
    city2ctry=np.array([ctry[city==i][0] for i in range(k*cities)])
    hits=0; ests=[]
    for _ in range(nsim):
        sd=np.where(M,np.sqrt(sh['mob'])[:,None],np.sqrt(sh['oth'])[:,None])  # 3 x G
        Y=(rng.normal(0,1,(k,NM+NO))*sd[0])[ctry]+(rng.normal(0,1,(k*cities,NM+NO))*sd[1])[city]+rng.normal(0,1,(n,NM+NO))*sd[2]+(rng.normal(0,np.sqrt(SEASON),(ns,NM+NO)))[season_key]
        o=r2_fast(Y,ctry,k); D=o[M].mean()-o[~M].mean(); ests.append(o.mean())
        cnt=0
        for _ in range(nperm):
            perm=rng.permutation(city2ctry)          # reassign whole cities to countries
            oo=r2_fast(Y,perm[city],k)
            if (oo[M].mean()-oo[~M].mean())>=D: cnt+=1
        if (cnt+1)/(nperm+1)<0.05: hits+=1
    return hits/nsim, float(np.mean(ests))
design=sys.argv[1]; out={}
for k in [8,12,16,24,32]:
    for c in [1,2,3,5]:
        p,est=simulate(k,c,design=design)
        out[f'{k}|{c}']=[round(p,2),round(est,3)]
        print(design,k,c,p,round(est,3),flush=True)
json.dump(out,open(f'design_{design}.json','w'),indent=1)
