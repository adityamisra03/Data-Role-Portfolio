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
# pip install matplotlib-venn

#######################
## Load up libraries ##
#######################

import pandas as pd
import GEOparse as geo
from scipy.stats import hypergeom
from statsmodels.stats import multitest
import gseapy as gp
import mygene
import os
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

##################
## Load up data ##
##################

CRISPR_raw = pd.read_csv('data_enrichment/CRISPR_TFs.csv', header=0, index_col=0)
CRISPR_131=[]
CRISPR_827=[]
for i in range(1445):
    if CRISPR_raw['0131_FDR'][i] <= 0.05:
        CRISPR_131.append(CRISPR_raw['Entrez ID'][i])
for i in range(1445):
    if CRISPR_raw['0827_FDR'][i] <= 0.05:
        CRISPR_827.append(CRISPR_raw['Entrez ID'][i])
CRISPR = set(CRISPR_827 + CRISPR_131)

DEG = (pd.read_csv('data_enrichment/DEG_TFs.csv',header=0, index_col=False))
DEG_E = set(DEG['Entrez ID'])

DepMap = (pd.read_csv('data_enrichment/DepMap_gbm.csv', header=None, index_col=False))
print(DepMap)
DepMap_E = set(set(CRISPR_raw['Entrez ID']).intersection(set(DepMap[0])))
print(DepMap_E)

Sygnal = (pd.read_csv('data_enrichment/gbmSYGNAL_TFs.csv', header=0, index_col=False))
SygnalTF = [j for i in Sygnal['Overlapping'] for j in i.split(' ')]
Entrez_Sygnal = set(CRISPR_raw.loc[SygnalTF, 'Entrez ID'])

gbmGenes = pd.read_csv('data_enrichment/DisGeNET_gbmGenes.txt', sep='\t', header=0, index_col=False) 
gbmGenes.dropna(axis=1, inplace=True)
gbmGenes_E = set(gbmGenes['EntrezGeneId'])

allTFs = list(CRISPR_raw['Entrez ID'])


#################################
## Hypergeometric comparisons: ##
## Primary Goal                ##
#################################

# Sygnal vs CRISPR (cutoff 0.01)
k = len(Entrez_Sygnal.intersection(CRISPR))
M = len(allTFs)
n = len(Entrez_Sygnal)
N = len(CRISPR)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

# Sygnal vs DepMap (cutoff 0.01)
k = len(Entrez_Sygnal.intersection(DepMap_E))
M = len(allTFs)
n = len(Entrez_Sygnal)
N = len(DepMap_E)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

# Sygnal vs DEG (cutoff 0.01)
k = len(Entrez_Sygnal.intersection(DEG_E))
M = len(allTFs)
n = len(Entrez_Sygnal)
N = len(DEG_E)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

# Sygnal vs DisGenNET (cutoff 0.01)
k = len(Entrez_Sygnal.intersection(gbmGenes_E))
M = len(allTFs)
n = len(Entrez_Sygnal)
N = len(gbmGenes_E)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])


#################################
## Hypergeometric comparisons: ##
## Secondary Goal              ##
#################################

# Toledo Crispr vs DisGeNET (cutoff 0.01)
k = len(CRISPR.intersection(gbmGenes_E))
M = len(allTFs)
n = len(CRISPR)
N = len(gbmGenes_E)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])
venn2([CRISPR, gbmGenes_E], ('Toledo CRISPR-Cas9 hit TFs', 'DisGeNET'))
plt.show()

# DepMap vs DisGeNET (cutoff 0.01)
k = len(DepMap_E.intersection(gbmGenes_E))
M = len(allTFs)
n = len(DepMap_E)
N = len(gbmGenes_E)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])
venn2([DepMap_E, gbmGenes_E], ('DepMap CRISPR-Cas9 hit TFs', 'DisGeNET'))
plt.show()

# Toledo Crispr vs DepMap (cutoff 0.01)
k = len(CRISPR.intersection(DepMap_E))
M = len(allTFs)
n = len(CRISPR)
N = len(DepMap_E)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])
venn2([CRISPR, DepMap_E], ('Toledo CRISPR-Cas9 hit TFs', 'DepMap CRISPR-Cas9 hit TFs'))
plt.show()


#################################
## Hypergeometric comparisons: ##
## Tertiary Goal               ##
#################################
venn3([CRISPR, DepMap_E, Entrez_Sygnal], ('Toledo CRISPR-Cas9 TF hits', 'DepMap CRISPR-Cas9 TF hits', 'gbmSYGNAL TFs'))
plt.show()

venn3([CRISPR, DepMap_E, gbmGenes_E], ('Toledo CRISPR-Cas9 TF hits', 'DepMap CRISPR-Cas9 TF hits', 'DisGeNET GBM TFs'))
plt.show()
