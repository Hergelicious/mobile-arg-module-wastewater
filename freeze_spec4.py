"""FROZEN SPECIFICATION 4 - clinical relevance. Written and hashed BEFORE any association is computed.

Question: does the wastewater module carry information about clinical antibiotic resistance, and does it add
anything beyond national antibiotic consumption and economic covariates?

Cohort: the external sewage cohort (module score available for 112 countries). The discovery cohort is not used:
only 8 of its 16 countries have GLASS data, and the USA and China are not among them.

Clinical outcome (WHO GLASS 2022 compilation, bloodstream isolates 2017-2020, pooled over years per country as
Resistant / InterpretableAST):
  primary   : Escherichia coli resistant to third-generation cephalosporins (cefotaxime, ceftriaxone or ceftazidime)
  secondary : E. coli resistant to fluoroquinolones (ciprofloxacin); Klebsiella pneumoniae resistant to
              third-generation cephalosporins; K. pneumoniae resistant to carbapenems (meropenem or imipenem)
Countries are included when at least 20 interpretable isolates are available for the phenotype.

Predictors (country means, sewage cohort): module score; acquired-but-not-module score; total ARG abundance; ARG richness.
Covariates: national antibiotic consumption (DDD per 1000 per day), log GDP per capita (PPP, WDI), WHO region.

Tests:
  T1 (primary)   Spearman correlation between module score and the primary outcome.
  T2 (adjustment) linear model logit(resistance) ~ module score + log consumption + log GDP per capita + region;
                  the question is whether the module coefficient remains different from zero.
  T3 (specificity) the same model with the acquired-but-not-module score, and with total ARG abundance, in place of
                  the module score; an association that is equally strong for those is not module-specific.
  T4 (added value) change in adjusted R2 when the module score is added to a model containing only consumption,
                  GDP and region.
Success: T1 P<0.05 in the positive direction AND the module coefficient in T2 P<0.05.
Uninformative if fewer than 25 countries have both wastewater and clinical data.
Every outcome is reported, including a null result, and T3 is reported whatever T1 and T2 show.
"""
import json, hashlib, re, numpy as np, pandas as pd
d=pd.read_csv('/tmp/glass/GLASS2022-master/compiled_WHO_GLASS_2022.csv')
meta=pd.read_excel('/mnt/user-data/uploads/41467_2025_66070_MOESM4_ESM.xlsx').drop_duplicates('complete_name')
nm=lambda s: re.sub(r'[^a-z]','',str(s).lower())
glass={nm(c) for c in d.CountryTerritoryArea.unique()}; sew={nm(c) for c in meta.country.dropna().unique()}
spec={'version':'frozen-4','cohort':'external sewage cohort (Martiny et al. 2025)',
 'clinical_source':'WHO GLASS 2022 compilation (compiled_WHO_GLASS_2022.csv), BLOOD specimens, 2017-2020',
 'outcomes':{'primary':'E. coli, third-generation cephalosporins (Cefotaxime/Ceftriaxone/Ceftazidime)',
             'secondary':['E. coli ciprofloxacin','K. pneumoniae third-generation cephalosporins','K. pneumoniae carbapenems (Meropenem/Imipenem)']},
 'pooling':'sum Resistant / sum InterpretableAST across years and listed drugs, per country; minimum 20 interpretable isolates',
 'predictors':['module score','acquired-not-module score','total ARG abundance','ARG richness'],
 'covariates':['antibiotic consumption (DDD/1000/day)','log GDP per capita PPP (WDI NY.GDP.PCAP.PP.KD, latest year <= 2021)','WHO region'],
 'tests':{'T1':'Spearman: module score vs primary outcome','T2':'logit(resistance) ~ module + log consumption + log GDP + region',
          'T3':'same model with acquired-not-module, and with total ARG, in place of module (specificity)',
          'T4':'change in adjusted R2 from adding the module to covariates only'},
 'success':'T1 P<0.05 positive AND module coefficient in T2 P<0.05','uninformative_if':'fewer than 25 countries with both data types',
 'countries_with_both_before_filtering':len(glass&sew),
 'reporting':'all outcomes reported, including null results; T3 reported regardless'}
txt=json.dumps(spec,indent=1,sort_keys=True); open('frozen_spec4.json','w').write(txt)
print('countries with wastewater and GLASS data:',len(glass&sew))
print('SHA-256:',hashlib.sha256(txt.encode()).hexdigest())
