## Functional enrichment analysis
#pip install mygene
#pip install matplotlib-venn

#######################
## Load up libraries ##
#######################

import pandas as pd
import GEOparse as geo
from scipy.stats import hypergeom
from statsmodels.stats import multitest
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

import mygene

##################
## Load up data ## want all transcription factors limit them by CRISPR_TFs
##################    list(set(list1+list2)) use

# Get total possible (background) gene set
CRISPR_TFs = pd.read_csv("data_enrichment/CRISPR_TFs.csv", header = 0, index_col = 0)
CRISPR_131 = []
CRISPR_827 = []
for i in range(1445):
    if CRISPR_TFs['0131_FDR'][i] <= 0.05:
        CRISPR_131.append(CRISPR_TFs['Entrez ID'][i])
for i in range(1445):
    if CRISPR_TFs['0827_FDR'][i] <= 0.05:
        CRISPR_827.append(CRISPR_TFs['Entrez ID'][i])
CRISPR = set(CRISPR_131 + CRISPR_827)
len(CRISPR) #435 (381 non-duplicates)
print(type(CRISPR))

DepMap = pd.read_csv("data_enrichment/DepMap_gbm.csv", header = None, index_col = False)
DepMap = set(set(CRISPR_TFs['Entrez ID']).intersection(set(DepMap[0])))
len(DepMap) #2338
print(type(DepMap))

Deg_TFs = pd.read_csv("data_enrichment/DEG_TFs.csv", header = 0, index_col=0)
Deg_TFs = set(Deg_TFs['Entrez ID'])
len(Deg_TFs) #264
print(type(Deg_TFs))

gbmSYGNAL_TFs = pd.read_csv("data_enrichment/gbmSYGNAL_TFs.csv", header = 0, index_col=0)
gbmSYGNAL_TFs = [j for i in gbmSYGNAL_TFs['Overlapping'] for j in i.split(' ')]
SYGNAL_TFs = set(CRISPR_TFs.loc[gbmSYGNAL_TFs, 'Entrez ID'])
len(SYGNAL_TFs) #74
print(type(SYGNAL_TFs))

DisGeNET_gbmGenes_all = pd.read_csv("data_enrichment/DisGeNET_gbmGenes.txt", sep = '\t', header = 0, index_col = False)
DisGeNET_gbmGenes_all.dropna(axis=1, inplace=True)
DisGeNET_gbmGenes = set(set(CRISPR_TFs['Entrez ID']).intersection(DisGeNET_gbmGenes_all['EntrezGeneId']))
len(DisGeNET_gbmGenes) #1247
print(type(DisGeNET_gbmGenes))

allTFs = list(CRISPR_TFs['Entrez ID'])
len(allTFs) #1445
print(type(allTFs))

#################################
## Hypergeometric comparisons: ##
## Primary Goal                ##
#################################

# gbmSYGNAL vs. Toledo CRISPR_TFs (cutoff 0.01)
k = len(SYGNAL_TFs.intersection(CRISPR))
M = len(allTFs)
n = len(SYGNAL_TFs)
N = len(CRISPR)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

# gbmSYGNAL vs. DepMap
k = len(SYGNAL_TFs.intersection(DepMap))
M = len(allTFs)
n = len(SYGNAL_TFs)
N = len(DepMap)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

# gbmSYGNAL vs. Deg_TFs
k = len(SYGNAL_TFs.intersection(Deg_TFs))
M = len(allTFs)
n = len(SYGNAL_TFs)
N = len(Deg_TFs)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

# gbmSYGNAL vs. DisGeNET_gbmGenes
k = len(SYGNAL_TFs.intersection(DisGeNET_gbmGenes))
M = len(allTFs)
n = len(SYGNAL_TFs)
N = len(DisGeNET_gbmGenes)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])

#################################
## Hypergeometric comparisons: ##
## Secondary Goal              ##
#################################

# Toledo Crispr vs DisGeNET (cutoff 0.01)
k = len(CRISPR.intersection(DisGeNET_gbmGenes))
M = len(allTFs)
n = len(CRISPR)
N = len(DisGeNET_gbmGenes)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])
venn2([CRISPR, DisGeNET_gbmGenes], ('Toledo CRISPR-Cas9 hit TFs', 'DisGeNET'))
plt.show()

# DepMap vs DisGeNET (cutoff 0.01)
k = len(DepMap.intersection(DisGeNET_gbmGenes))
M = len(allTFs)
n = len(DepMap)
N = len(DisGeNET_gbmGenes)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])
venn2([DepMap, DisGeNET_gbmGenes], ('DepMap CRISPR-Cas9 hit TFs', 'DisGeNET'))
plt.show()

# Toledo Crispr vs DepMap (cutoff 0.01)
k = len(CRISPR.intersection(DepMap))
M = len(allTFs)
n = len(CRISPR)
N = len(DepMap)
p_value = hypergeom.sf(k, M, n, N)
print([k,M,n,N,p_value])
venn2([CRISPR, DepMap], ('Toledo CRISPR-Cas9 hit TFs', 'DepMap CRISPR-Cas9 hit TFs'))
plt.show()

# Entrez Id - Gene symbol conversion

CRISPR_TFs_geneConv = pd.read_csv("data_enrichment/CRISPR_TFs.csv", header = 0, index_col = False)
mapping = dict(zip(CRISPR_TFs_geneConv['Entrez ID'], CRISPR_TFs_geneConv['TF']))

DisGeNET_gbmGenes_all = DisGeNET_gbmGenes_all.rename(columns={"EntrezGeneId": "Entrez ID"})
DisGeNET_gbmGenes_df = pd.DataFrame(DisGeNET_gbmGenes, columns=['Entrez ID'])
DisGeNET_gbmGenes_gene = set(DisGeNET_gbmGenes_df['Entrez ID'].map(mapping))

CRISPR_df = pd.DataFrame(CRISPR, columns=['Entrez ID'])
CRISPR_TFs_gene = set(CRISPR_df['Entrez ID'].map(mapping))

DepMap_gene_df = pd.DataFrame(DepMap, columns=['Entrez ID'])
DepMap_gene = set(DepMap_gene_df['Entrez ID'].map(mapping))

gbmSYGNAL_TFs = set(gbmSYGNAL_TFs)

# make sure compatible for venn non-gene conversion set
CRISPR = set(CRISPR)
DepMap = set(DepMap)
SYGNAL_TFs = set(SYGNAL_TFs)
DisGeNET_gbmGenes = set(DisGeNET_gbmGenes)

#################################
## Hypergeometric comparisons: ##
## Tertiary Goal               ##
#################################

# Comparing three sets Entrez IDs so totals
v = venn3([CRISPR, DepMap, SYGNAL_TFs], ('Toledo CRISPR-Cas9 TF hits', 'DepMap CRISPR-Cas9 TF hits', 'gbmSYGNAL TFs'))
plt.show()


# Comparing three sets Entrez IDs so totals
venn3([CRISPR, DepMap, DisGeNET_gbmGenes], ('Toledo CRISPR-Cas9 TF hits', 'DepMap CRISPR-Cas9 TF hits', 'DisGeNET GBM TFs'))
plt.show()



# Comparing three sets gene conversion included only including known genes from the CRISPR data so much smaller numbers - lots of nan for genes
venn3([CRISPR_TFs_gene, DepMap_gene, gbmSYGNAL_TFs], ('Toledo CRISPR-Cas9 TF hits', 'DepMap CRISPR-Cas9 TF hits', 'gbmSYGNAL TFs'))
plt.show()

print("Overlap between CRISPR and DepMap:", CRISPR_TFs_gene.intersection(DepMap_gene))
CRISPR_DepMap = pd.DataFrame(CRISPR_TFs_gene.intersection(DepMap_gene), columns=['Overlap between CRISPR and DepMap'])

print("Overlap between CRISPR and SYGNAL_TFs:", CRISPR_TFs_gene.intersection(gbmSYGNAL_TFs))
CRISPR_SYGNAL_TFs = pd.DataFrame(CRISPR_TFs_gene.intersection(gbmSYGNAL_TFs), columns=['Overlap between CRISPR and SYGNAL_TFs'])

print("Overlap between DepMap and SYGNAL_TFs:", DepMap_gene.intersection(gbmSYGNAL_TFs))
DepMap_SYGNAL_TFs = pd.DataFrame(DepMap_gene.intersection(gbmSYGNAL_TFs), columns=['Overlap between DepMap and SYGNAL_TFs'])

print("Overlap between all three sets:", CRISPR_TFs_gene.intersection(DepMap_gene, gbmSYGNAL_TFs))
CRISPR_DepMap_SYGNAL_TFs = pd.DataFrame(CRISPR_TFs_gene.intersection(DepMap_gene, gbmSYGNAL_TFs), columns=['Overlap between all three sets'])

CRISPR_DepMap_SYGNAL_Overlap_Genes = pd.concat([CRISPR_DepMap, CRISPR_SYGNAL_TFs, DepMap_SYGNAL_TFs, CRISPR_DepMap_SYGNAL_TFs])
CRISPR_DepMap_SYGNAL_Overlap_Genes.to_csv('CRISPR_DepMap_SYGNAL_Overlap_Genes.csv')

# Comparing three sets sets gene conversion included only including known genes from the CRISPR data so much smaller numbers - lots of nan for genes
venn3([CRISPR_TFs_gene, DepMap_gene, DisGeNET_gbmGenes_gene], ('Toledo CRISPR-Cas9 TF hits', 'DepMap CRISPR-Cas9 TF hits', 'DisGeNET GBM TFs'))
plt.show()

print("Overlap between CRISPR and DepMap:", CRISPR_TFs_gene.intersection(DepMap_gene))

print("Overlap between CRISPR and DisGeNET_gbmGenes_gene:", CRISPR_TFs_gene.intersection(DisGeNET_gbmGenes_gene))
CRISPR_DisGeNET_gbmGenes = pd.DataFrame(CRISPR_TFs_gene.intersection(DisGeNET_gbmGenes_gene), columns=['Overlap between CRISPR and DisGeNET_gbmGenes'])

print("Overlap between DepMap and DisGeNET_gbmGenes_gene:", DepMap_gene.intersection(DisGeNET_gbmGenes_gene))
DepMap_DisGeNET_gbmGenes = pd.DataFrame(DepMap_gene.intersection(DisGeNET_gbmGenes_gene), columns=['Overlap between DepMap and DisGeNET_gbmGenes'])

print("Overlap between all three sets:", CRISPR_TFs_gene.intersection(DepMap_gene, DisGeNET_gbmGenes_gene))
CRISPR_DepMap_DisGeNET_gbmGenes = pd.DataFrame(CRISPR_TFs_gene.intersection(DepMap_gene, DisGeNET_gbmGenes_gene), columns=['Overlap between all three sets'])

CRISPR_DepMap_DisGeNET_gbmGenes_Overlap_Genes = pd.concat([CRISPR_DepMap, CRISPR_DisGeNET_gbmGenes, DepMap_DisGeNET_gbmGenes, CRISPR_DepMap_DisGeNET_gbmGenes])
CRISPR_DepMap_DisGeNET_gbmGenes_Overlap_Genes.to_csv('CRISPR_DepMap_DisGeNET_gbmGenes_Overlap_Genes.csv')