"""(A) Frozen surveillance panel: the module and mobile-associated gene sets with cross-database names."""
import json, re, pandas as pd
g=pd.read_csv('gene_table_final.csv',index_col=0)
mm=pd.read_csv('module_membership.csv').set_index('gene').loc[g.index]
S3=json.load(open('frozen_spec3.json')); S2=json.load(open('frozen_spec2.json'))
mapped=S3['gene_mapping']['module_mapped']
def norm_card(x):
    y=x.strip()
    m=re.fullmatch(r'(Erm|Ere|Mph|Msr|Lnu|Mef|Lsa|Vat|Vga)([A-Z])(\d*)',y,re.I)
    if m: return f"{m.group(1).lower()}({m.group(2).lower()}){m.group(3)}"
    m=re.fullmatch(r'(OXA|TEM|GES|VEB|PER|CARB|CTX-M|SHV|IMP|VIM|NDM|KPC|CMY|ACT|LCR|NPS)-?(\d+)',y,re.I)
    if m: return f"bla{m.group(1).lower()}-{m.group(2)}"
    return y.lower()
rf=set(S2['D_external_classification']['resfinder_acquired_genes'])
p=pd.DataFrame({'gene_CARD':g.index,
                'in_module':mm.in_module.values,
                'mobile_associated_curated':g.M_prior.values,
                'resfinder_acquired':[x in rf for x in g.index],
                'name_normalised':[norm_card(x) for x in g.index],
                'panres_ids':[';'.join(mapped.get(x,[])) for x in g.index],
                'drug_class':g['class'].values,'mechanism':g.mechanism.values,
                'prevalence_discovery':g.prevalence.round(3).values,
                'country_R2_discovery':g.R2_country.round(3).values,
                'rho_intI1_within_country':pd.read_csv('intI1_rho.csv',index_col=0).loc[g.index].iloc[:,0].round(3).values})
p.to_csv('resistome_module_panel.csv',index=False)
print('panel rows',len(p),'| module',int(p.in_module.sum()),'| mobile-associated',int(p.mobile_associated_curated.sum()),'| with PanRes id',int((p.panres_ids!='').sum()))
