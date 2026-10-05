"""FROZEN SPECIFICATION 5 - within-country (panel) and drug-class-matched tests of the clinical link.
Written and hashed BEFORE any association is computed. Both tests address limitations stated in the manuscript:
the clinical comparison so far is cross-sectional (confounded by national income) and uses the module as a whole
(not matched to the phenotype it should mechanistically predict).

TEST P (within-country panel). Country-year sewage measures against country-year clinical resistance, with country
and year fixed effects, so every stable national difference (income, health system, laboratory practice) is removed.
  outcome    : logit of resistant/interpretable per country-year (>=20 interpretable isolates), E. coli against
               third-generation cephalosporins (primary); the four phenotypes pooled with phenotype fixed effects
               (secondary, to gain observations)
  predictors : module score; acquired-but-not-module score; total ARG abundance (country-year means, ALR)
  model      : OLS with country and year fixed effects, cluster-robust standard errors by country
  POWER      : only 11 countries have two or more matched years (50 country-years). This test is therefore
               declared LOW-POWERED IN ADVANCE: a null result will be reported as inconclusive, not as evidence of
               absence, and the estimate will be reported with its confidence interval.

TEST C (drug-class matched). Instead of the whole module, the ARG set is matched to the phenotype:
  E. coli third-generation cephalosporins  <- beta_lactam ARGs
  E. coli ciprofloxacin                    <- quinolone ARGs
  K. pneumoniae third-generation cephalosporins <- beta_lactam ARGs
  K. pneumoniae carbapenems                <- carbapenemase genes (bla NDM, KPC, OXA-48-like, VIM, IMP)
  statistic : Spearman correlation across countries, and the coefficient in a model adjusting for log GDP per
              capita and WHO region; module-matched genes are additionally compared with the same class restricted
              to module members, to ask whether the module adds anything within the matched class
  multiplicity: four pairings, Holm correction on the primary Spearman tests
  success   : Holm-adjusted P<0.05 in the positive direction for at least one pairing
Both tests are reported whatever the outcome.
"""
import json, hashlib
spec={'version':'frozen-5',
 'test_P':{'design':'country and year fixed effects, cluster-robust SE by country','outcome':'logit(resistant/interpretable) per country-year, >=20 isolates',
   'primary_outcome':'E. coli third-generation cephalosporins','secondary':'four phenotypes pooled with phenotype fixed effects',
   'predictors':['module score','acquired-not-module score','total ARG abundance'],
   'power':'11 countries with >=2 matched years, 50 country-years; declared low-powered in advance; a null is inconclusive'},
 'test_C':{'pairings':{'ecoli_3GC':'beta_lactam','ecoli_cipro':'quinolone','kpneu_3GC':'beta_lactam','kpneu_carb':'carbapenemase (NDM/KPC/OXA-48-like/VIM/IMP)'},
   'statistics':['Spearman across countries','coefficient adjusting for log GDP per capita and WHO region','module-restricted version of the matched class'],
   'multiplicity':'Holm across the four primary Spearman tests','success':'Holm-adjusted P<0.05 positive for at least one pairing'},
 'reporting':'both tests reported whatever the outcome; no test is redefined after the data are seen'}
txt=json.dumps(spec,indent=1,sort_keys=True); open('frozen_spec5.json','w').write(txt)
print('SHA-256:',hashlib.sha256(txt.encode()).hexdigest())
