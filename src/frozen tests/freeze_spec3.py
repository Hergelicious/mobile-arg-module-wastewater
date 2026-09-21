"""FROZEN SPECIFICATION 3 - Tier A tests in external data. Written and hashed BEFORE the data files are obtained.
Gene mapping is built mechanically from Martiny et al. Supplementary Data 2 (already in hand); no outcome values exist yet."""
import json, re, hashlib, pandas as pd
g=pd.read_csv('gene_table_final.csv',index_col=0); core=list(g.index)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[core]
MOD=[x for x in core if bool(mm.loc[x,'in_module'])]
sd2=pd.read_excel('/mnt/user-data/uploads/41467_2025_66070_MOESM5_ESM.xlsx')
def norm_card(x):                      # identical to frozen specs 1 and 2
    y=x.strip()
    m=re.fullmatch(r'(Erm|Ere|Mph|Msr|Lnu|Mef|Lsa|Vat|Vga)([A-Z])(\d*)',y,re.I)
    if m: return f"{m.group(1).lower()}({m.group(2).lower()}){m.group(3)}"
    m=re.fullmatch(r'(OXA|TEM|GES|VEB|PER|CARB|CTX-M|SHV|IMP|VIM|NDM|KPC|CMY|ACT|LCR|NPS)-?(\d+)',y,re.I)
    if m: return f"bla{m.group(1).lower()}-{m.group(2)}"
    return y.lower()
acq=sd2[sd2.database=='resfinder'].copy()
acq['stem']=acq.fa_name.astype(str).str.split('|',n=1).str[1].str.replace(r'_\d+_[^_]+$','',regex=True).str.lower()
stem2pan=acq.groupby('stem').gene.apply(lambda s: sorted(set(s))).to_dict()
mapped={x:stem2pan[norm_card(x)] for x in core if norm_card(x) in stem2pan}
mod_mapped={x:v for x,v in mapped.items() if x in MOD}
fg=sorted(set(sd2[sd2.database!='resfinder'].gene))
spec={'version':'frozen-3','purpose':'Tier A: external validation, temporal stability, drivers and cross-compartment persistence',
 'inputs_expected':['panres_counts.csv and kingdom_and_motus_agg_counts.csv (Zenodo 14652833)',
                    'SraRunTable for PRJNA509305 (GWMC run/batch metadata)',
                    'Browne et al. 2021 national antibiotic consumption (DDD per 1000 per day)'],
 'normalisation':'length-adjusted ARG fragment counts, additive log-ratio with total bacterial (mOTU) fragments as the reference component, as in Martiny et al.; zeros imputed as 65% of the per-sample detection limit; city means used as the unit for geographic tests',
 'gene_mapping':{'rule':'CARD name normalised by the shared function, matched to the stem of a ResFinder entry in Supplementary Data 2, then to its PanRes gene id(s)',
   'n_core_mapped':len(mapped),'n_module_mapped':len(mod_mapped),'module_mapped':mod_mapped,
   'comparison_acquired_not_module':sorted(set(acq.gene)-{p for v in mod_mapped.values() for p in v}),
   'n_functional_genes':len(fg)},
 'tests':{
  'D_validation':{'primary':'D_ext = mean marginal country R2 of mapped module genes minus mean of acquired genes not in the module, city means, cities nested in countries',
    'null':'9,999 permutations reassigning whole cities among countries','success':'P<0.05, positive direction',
    'secondary':['same statistic with functional (latent) genes as the comparison set - expected positive and largely a replication of Martiny et al.',
                 'within-country mean pairwise correlation of mapped module genes vs 1,000 abundance-matched random sets of other acquired genes'],
    'uninformative_if':'fewer than 15 module genes map'},
  'E_temporal':{'question':'is the country signal stable across sampling years, i.e. not an artefact of one campaign',
    'primary':'Spearman correlation across cities sampled in >=2 different years between their module score in the earliest and latest year; success P<0.05',
    'secondary':['variance components of the module score for country, city, year and residual (REML); success if the country component exceeds the year component',
                 'D_ext recomputed within each sampling year separately'],
    'uninformative_if':'fewer than 30 cities have two or more years'},
  'F_drivers':{'primary':'Spearman correlation across countries between the mean module score and national antibiotic consumption (DDD/1000/day, nearest available year, capped at 2018)',
    'secondary':['partial correlation controlling for log GDP per capita','mixed model with world region as a random intercept'],
    'success':'P<0.05, positive direction','note':'Njage et al. 2023 modelled socio-economic drivers of total sewage AMR; this test is specific to the module'},
  'G_cross_compartment':{'question':'does a city sewage signature predict its activated-sludge signature',
    'primary':'Spearman correlation across the 12 cities present in both cohorts between the module share in Martiny sewage and in GWMC activated sludge',
    'success':'P<0.05','status':'exploratory: n=12, different years, different pipelines; reported regardless'},
  'H_batch':{'question':'does the GWMC country effect survive adjustment for sequencing run and submission batch',
    'primary':'D (mobile-associated vs other, GWMC) recomputed after residualising each gene on run-level factors extracted from the SRA run table',
    'success':'D remains positive with P<0.05','note':'campaign dates are absent from the published metadata; only run-level fields can be used'}},
 'reporting':'every test above is reported in the manuscript whatever its outcome; no test is redefined after the data are seen'}
txt=json.dumps(spec,indent=1,sort_keys=True); open('frozen_spec3.json','w').write(txt)
print('core genes mapped to PanRes:',len(mapped),'| module genes mapped:',len(mod_mapped))
print('SHA-256:',hashlib.sha256(txt.encode()).hexdigest())
