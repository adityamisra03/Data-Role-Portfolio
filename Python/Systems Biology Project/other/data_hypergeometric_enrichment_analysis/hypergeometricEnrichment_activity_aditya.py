##########################################################
## OncoMerge:  hypergeometricEnrichment.py              ##
##  ______     ______     __  __                        ##
## /\  __ \   /\  ___\   /\ \/\ \                       ##
## \ \  __ \  \ \___  \  \ \ \_\ \                      ##
##  \ \_\ \_\  \/\_____\  \ \_____\                     ##
##   \/_/\/_/   \/_____/   \/_____/                     ##
## @Developed by: Plaisier Lab                          ##
##   (https://plaisierlab.engineering.asu.edu/)         ##
##   Arizona State University                           ##
##   242 ISTB1, 550 E Orange St                         ##
##   Tempe, AZ  85281                                   ##
## @Author:  Chris Plaisier                             ##
## @License:  GNU GPLv3                                 ##
##                                                      ##
## If this program is used in your analysis please      ##
## mention who built it. Thanks. :-)                    ##
##########################################################

#########################
## Will need to install ##
#########################

## Functional enrichment analysis
# pip install mygene
# pip install gseapy

#######################
## Load up libraries ##
#######################

import pandas as pd
import GEOparse as geo
from scipy.stats import hypergeom
from statsmodels.stats import multitest
import gseapy as gp
import mygene
mg = mygene.MyGeneInfo()

##################
## Load up data ##
##################

# Get total possible (background) gene set
df_CRISPR_TFs = pd.read_csv('data_enrichment/CRISPR_TFs.csv', header=0, index_col=0)
df_DEG_TFs = pd.read_csv('data_enrichment/DEG_TFs.csv', header=0, index_col=False)
df_DepMap_gbm = pd.read_csv('data_enrichment/DepMap_gbm.csv', header=None, index_col=False)
df_gbmSYGNAL_TFs = pd.read_csv('data_enrichment/gbmSYGNAL_TFs.csv', header=0, index_col=False)
Sygnal = [j for i in df_gbmSYGNAL_TFs['Overlapping'] for j in i.split(' ')]
gbmGenes = pd.read_csv('data_enrichment/DisGeNet_gbmGenes.txt', sep='\t', header=0, index_col=False)
allTFs = list(df_CRISPR_TFs['Entrez ID'])
Entrez = df_CRISPR_TFs.loc(Sygnal, 'Entrez ID')


CRISPR_TFs_0131_FDR = set(df_CRISPR_TFs['0131_FDR'])
CRISPR_TFs_0827_FDR = set(df_CRISPR_TFs['0827_FDR'])
filt_CRISPR_TFs_0131_FDR = set([i for i in CRISPR_TFs_0131_FDR if i <= 0.05])
filt_CRISPR_TFs_0827_FDR = set([i for i in CRISPR_TFs_0827_FDR if i <= 0.05])
print(filt_CRISPR_TFs_0131_FDR)
print(filt_CRISPR_TFs_0827_FDR)

# There are some funky IDs with ///, so lets clean that up by taking splitting them and keeping all IDs
# Uses very helpful list comprehension:  https://realpython.com/list-comprehension-python/
# Also uses a helpful trick in list comprehension for un-nesting lists:  [item for sublist in l for item in sublist]
allGenes = set([int(i3) for i2 in [i.split(' /// ') for i in allGenes] for i3 in i2])
print(len(allGenes)) # 9,946

# Load up USF1 over-expression (OE) differentially expressed genes (DEGs)
# https://doi.org/10.1371/journal.pgen.1000642.s003
df_USF1_OE_DEGs = pd.read_csv('data_enrichment/df_USF1_OE_DEGs.csv', header=0, index_col=0) 
USF1_OE_DEGs = set(df_USF1_OE_DEGs['Entrez ID'].dropna()).intersection(allGenes)
print(len(USF1_OE_DEGs)) # 2140

# Load up USF1 target genes (from TFBS_DB)
# http://tfbsdb.systemsbiology.net/searchtf?searchterm=V_USF_02_M00122
# http://tfbsdb.systemsbiology.net/tfgenes_csv/V_USF_02_M00122
df_USF_M00122_target_genes = pd.read_csv('data_enrichment/V_USF_02_M00122_genes.tsv', delimiter='\t', index_col=0, header=0)
USF_M00122_target_genes = set([i for i in df_USF_M00122_target_genes.index]).intersection(allGenes)
print(len(USF_M00122_target_genes)) # 319 (NEW source)

# Load up rs3737787 correlated genes
# https://doi.org/10.1371/journal.pgen.1000642.s004
df_rs3737787_corr = pd.read_csv('data_enrichment/df_rs3737787_corr.csv', header=0, index_col=0)
df_rs3737787_corr_up = df_rs3737787_corr.loc[df_rs3737787_corr['Effect \nEstimate'].astype(float)>0]
df_rs3737787_corr_up_0_01 = df_rs3737787_corr_up.loc[df_rs3737787_corr_up['P-value'].astype(float)<0.01]
rs3737787_corr_up_0_01 = set(df_rs3737787_corr_up_0_01['Entrez ID'].dropna()).intersection(allGenes)
print(len(rs3737787_corr_up_0_01)) # 145 (very different than what is in the paper, because we are using extra filters)

# Load up USF1 over-expression (OE) differentially expressed genes (DEGs)
# https://doi.org/10.1371/journal.pgen.1000642.s005
df_FCHL_DEGs = pd.read_csv('data_enrichment/df_FCHL_DEGs.csv', header=0, index_col=0)
FCHL_DEGs = set(df_FCHL_DEGs['Entrez ID'].dropna()).intersection(allGenes)
print(len(FCHL_DEGs)) # 2058

#################################
## Hypergeometric comparisons: ##
## rs3737787 is linked to USF1 ##
#################################

# USF1 OE DEGs vs. rs3737787 correlated/upregulated genes (cutoff 0.01)
k = len(USF1_OE_DEGs.intersection(rs3737787_corr_up_0_01))
M = len(allGenes)
n = len(USF1_OE_DEGs)
N = len(rs3737787_corr_up_0_01)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

# USF1 OE DEGs vs. USF1 target genes
k = len(USF_M00122_target_genes.intersection(USF1_OE_DEGs))
M = len(allGenes)
n = len(USF_M00122_target_genes)
N = len(USF1_OE_DEGs)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

# USF1 target genes vs. rs3737787 correlated genes
k = len(USF_M00122_target_genes.intersection(rs3737787_corr_up_0_01))
M = len(allGenes)
n = len(USF_M00122_target_genes)
N = len(rs3737787_corr_up_0_01)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

# FCHL DEGs vs. rs3737787 correlated genes
k = len(FCHL_DEGs.intersection(rs3737787_corr_up_0_01))
M = len(allGenes)
n = len(FCHL_DEGs)
N = len(rs3737787_corr_up_0_01)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])


#################################
## Hypergeometric comparisons: ##
## Clusters vs. gene-lists     ##
#################################

## Load up module memberships
module_membership = pd.read_csv('data_enrichment/moduleMembership_Plaisier2009.csv',header=0, index_col=0)
# Dictionary comprehension is like a for loop but the result is a key:value pair,
# useful for making a dictionary of modules to probe IDs
# https://towardsdatascience.com/10-examples-to-master-python-dictionary-comprehensions-7aaa536f5960
modules = {i:list(module_membership.loc[module_membership['WGCNA Module']==i].index) for i in set(module_membership['WGCNA Module'])}
print(modules.keys())
print(modules['Tan'])

# Now we need to convert the probe IDs into entrez IDs for the comparisons
# We will use the GPL570 from the fchl GSE we loaded earlier
gpl570 = fchl.gpls['GPL570'].table
gpl570.index = gpl570['ID']
# Nesting a list comprehension in a dictionary comprehension
modules_entrez = {i:set([int(i3) for i2 in [i.split(' /// ') for i in list(gpl570.loc[modules[i],'ENTREZ_GENE_ID'].dropna())] for i3 in i2]) for i in modules}
print(modules_entrez['Tan'])


## Compare all modules vs. FCHL DEGs
module_vs_FCHL = pd.DataFrame(index=modules_entrez.keys(),columns=['k','M','n','N','P-value'])
for module in modules_entrez:
    k = len(FCHL_DEGs.intersection(modules_entrez[module]))
    M = len(allGenes)
    n = len(FCHL_DEGs)
    N = len(modules_entrez[module])
    p_value = hypergeom.sf(k, M, n, N)
    module_vs_FCHL.loc[module] = [k,M,n,N,p_value]

# Add multiply hypothesis correction
module_vs_FCHL['Adj. P-value'] = multitest.multipletests(list(module_vs_FCHL['P-value']))[1]

# Print result
print(module_vs_FCHL.sort_values('Adj. P-value'))

# Write out CSV file
module_vs_FCHL.sort_values('Adj. P-value').to_csv('modules_vs_FCHL_DEGs.csv')


## Compare all modules vs. rs3737787 correlated/upregulated genes (cutoff 0.01)
module_vs_rs3737787 = pd.DataFrame(index=modules_entrez.keys(),columns=['k','M','n','N','P-value'])
for module in modules_entrez:
    k = len(rs3737787_corr_up_0_01.intersection(modules_entrez[module]))
    M = len(allGenes)
    n = len(rs3737787_corr_up_0_01)
    N = len(modules_entrez[module])
    p_value = hypergeom.sf(k, M, n, N)
    module_vs_rs3737787.loc[module] = [k,M,n,N,p_value]

# Add multiply hypothesis correction
module_vs_rs3737787['Adj. P-value'] = multitest.multipletests(list(module_vs_rs3737787['P-value']))[1]

# Print result
print(module_vs_rs3737787.sort_values('Adj. P-value'))

# Write out CSV file
module_vs_rs3737787.sort_values('Adj. P-value').to_csv('modules_vs_rs3737787_corr_up_0_01.csv')


#####################################
## Functional enrichment analysis: ##
## Characterizing the tan module   ##
#####################################

# Map entrez gene IDs to gene symbols
tan_mapped = mg.querymany(list(modules_entrez['Tan']), scopes='entrezgene', species=9606, as_dataframe=True)
tan_symbols = list(tan_mapped['symbol'].dropna())
background_mapped = mg.querymany(list(allGenes), scopes='entrezgene', species=9606, as_dataframe=True)
background_symbols = list(background_mapped['symbol'].dropna())

# Possible curated gene association libraries
names = gp.get_library_name()
print(names)

# Functional enrichment
enr1 = gp.enrichr(gene_list=tan_symbols, gene_sets=['GO_Biological_Process_2023','KEGG_2021_Human'], organism='Human', description='Tan_module', outdir='Tan_functional_enrichment', cutoff=0.05, no_plot=True) #, background=background_symbols) # tan_symbols
print(enr1.results)
# Check the Tan_functional_enrichment directory

