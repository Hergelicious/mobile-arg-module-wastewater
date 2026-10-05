"""Upgraded design simulation: more replicates and a seasonal-variance grid."""
import warnings; warnings.filterwarnings("ignore")
import json, sys, numpy as np
rng=np.random.default_rng(20260926)
K=json.load(open('canonical.json'))
sh={'mob':np.array([K['SS_country'][0],K['SS_city'][0],1-K['SS_country'][0]-K['SS_city'][0]]),
    'oth':np.array([K['SS_country'][1],K['SS_city'][1],1-K['SS_country'][1]-K['SS_city'][1]])}
NM,NO=45,69; M=np.array([True]*NM+[False]*NO)
def r2_fast(Y,lab,ng):
    cnt=np.bincount(lab,minlength=ng).astype(float); s=np.zeros((ng,Y.shape[1])); np.add.at(s,lab,Y)
    mu=s/cnt[:,None]; Yc=Y-Y.mean(0); return 1-((Y-mu[lab])**2).sum(0)/(Yc**2).sum(0)
def simulate(k,cities,season,design,plants=2,nsim=200,nperm=99):
    ctry=np.repeat(np.arange(k),cities*plants); city=np.repeat(np.arange(k*cities),plants); n=len(ctry)
    skey=ctry if design=='aliased' else np.tile(np.arange(4),int(np.ceil(n/4)))[:n]; ns=skey.max()+1
    c2c=np.array([ctry[city==i][0] for i in range(k*cities)])
    sd=np.where(M,np.sqrt(sh['mob'])[:,None],np.sqrt(sh['oth'])[:,None]); hits=0; ests=[]
    for _ in range(nsim):
        Y=(rng.normal(0,1,(k,NM+NO))*sd[0])[ctry]+(rng.normal(0,1,(k*cities,NM+NO))*sd[1])[city]+rng.normal(0,1,(n,NM+NO))*sd[2]+(rng.normal(0,np.sqrt(season),(ns,NM+NO)))[skey]
        o=r2_fast(Y,ctry,k); D=o[M].mean()-o[~M].mean(); ests.append(o.mean()); cnt=0
        for _ in range(nperm):
            oo=r2_fast(Y,rng.permutation(c2c)[city],k)
            if (oo[M].mean()-oo[~M].mean())>=D: cnt+=1
        if (cnt+1)/(nperm+1)<0.05: hits+=1
    return hits/nsim,float(np.mean(ests))
mode=sys.argv[1]; out={}
S0=json.load(open('frozen_results3_alr.json'))['E_temporal']['variance_fractions']['year']
if mode=='grid':
    for design in ['aliased','crossed']:
        for k in [8,12,16,24,32]:
            for c in [1,2,3,5]:
                p,e=simulate(k,c,S0,design); out[f'{design}|{k}|{c}']=[round(p,3),round(e,3)]
                print(design,k,c,p,flush=True)
else:
    for season in [0.15,0.25,0.35,0.45]:
        for c in [1,2,3,5]:
            p,e=simulate(16,c,season,'crossed'); out[f'season{season}|{c}']=[round(p,3),round(e,3)]
            print('season',season,c,p,flush=True)
json.dump(out,open(f'design2_{mode}.json','w'),indent=1)
