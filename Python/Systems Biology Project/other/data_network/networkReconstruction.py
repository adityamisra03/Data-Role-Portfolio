##########################################################
## OncoMerge:  networkReconstruction.py                 ##
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

#####################
## Need to install ##
#####################

## Install pyvis to export dynamic html visualizations of networks
# pip install pyvis
# pip install networkx


#######################
## Load up libraries ##
#######################

import pandas as pd
from pyvis.network import Network
import networkx as nx

##################
## Load up data ##
##################

# Load up SuppTab8 which has all the information needed
st8 = pd.read_csv('data_network/gbmSYGNAL.csv', header=0, index_col=0)

# Split out cancer from cluster name
st8['cancer'] = [i.split(' ')[1] for i in st8.index]

# Split up GO terms into a list
st8['GO Terms'] = [i.split('|') for i in st8['GO Terms']]

# Split up miRNAs into a list
st8['Mature Seqeunce IDs'] = [i.split('|') for i in st8['Mature Seqeunce IDs']]


## Load up accessory data
# Load up mature sequence IDs to miRNA names
miRDb = pd.read_csv('data_network/hsa_mature.csv', header=0, index_col=1)

# Convert GO terms into names
#gene2go = pd.read_csv('data/gene2go.hsa', delimiter='\t', header=None)
#go2name = gene2go.loc[gene2go[7]=='Process',[2,5]].drop_duplicates()
#go2name.index = go2name[2]


#################################################
## Filter network into cross-cancer subnetwork ##
#################################################

## Association to hallmark of cancer
hallmarks = ['Self Sufficiency In Growth Signals', 'Insensitivity To Antigrowth Signals', 'Evading Apoptosis', 'Limitless Replicative Potential', 'Sustained Angiogenesis', 'Tissue Invasion And Metastasis', 'Genome Instability And Mutation', 'Tumor Promoting Inflammation', 'Reprogramming Energy Metabolism', 'Evading Immune Detection']

# Function to select out significant hallmarks
# Takes as input a pandas series of hallmark similarity scores
def getHallmarks(scores, cutoff=0.8):
    return list((scores[scores>cutoff]).index)

# Run through all clusters and compute hallmarks
clusterHallmarks = []
for cluster in st8.index:
    scores = st8.loc[cluster, hallmarks]
    clusterHallmarks.append(getHallmarks(scores))

# Add to st8 dataFrame
st8['sigHallmarks'] = clusterHallmarks

# Subset st8 dataFrame to those with significant hallmarks
st8_hm = st8.loc[[True if len(st8.loc[cluster,'sigHallmarks'])>0 else False for cluster in st8.index]]


## Same oncogenic process regulated by same miRNA across cancers
# Occurences of miRNA in cancers
miRNAs = {}
for cluster in st8_hm.index:
    #print(row1)
    for miRNA in st8_hm.loc[cluster,'Mature Seqeunce IDs']:
        if not miRNA in miRNAs:
            miRNAs[miRNA] = []
        miRNAs[miRNA].append(cluster)

## Find cross-cancer miRNAs
cc_miRNAs = {}
cc_go_miRNAs = {}
for miRNA in miRNAs:
    # If have more than one cluster regulated by the miRNA
    if len(miRNAs[miRNA])>1:
        # If more than one cancer is regulated by the miRNA
        if len(st8_hm.loc[miRNAs[miRNA],'cancer'])>1:
            # Idenitfy which clusters have same go term -> hallmark relationships
            hm_dict = {}
            go_dict = {}
            # First, fill hm_dict with hallmarks and go_dict with go terms
            for cluster in miRNAs[miRNA]:
                # Build clusters so can test for common hallmarks
                for hm1 in st8_hm.loc[cluster,'sigHallmarks']:
                    if not hm1 in hm_dict:
                        hm_dict[hm1] = []
                    hm_dict[hm1].append(cluster)
                # Build clusters so can test for common GO terms
                for go1 in st8_hm.loc[cluster,'GO Terms']:
                    if not go1 in go_dict:
                        go_dict[go1] = []
                    go_dict[go1].append(cluster)
            # Second, find hallmarks in common across clusters, and refine subnetwork (cc_miRNAs)
            for hm1 in hm_dict:
                # If more than one cluster linked to same hallmark
                if len(hm_dict[hm1])>1:
                    if not miRNA in cc_miRNAs:
                        cc_miRNAs[miRNA] = {}
                    cc_miRNAs[miRNA][hm1] = miRNAs[miRNA]
            # Thrid, find go terms in common across clusters, and refine subnetwork
            for go1 in go_dict:
                # If more than one cluster linked to same GO term
                if len(go_dict[go1])>1:
                    # If had more than one cluster linked to same hallmark
                    if miRNA in cc_miRNAs:
                        if not miRNA in cc_go_miRNAs:
                            cc_go_miRNAs[miRNA] = {'hallmarks':cc_miRNAs[miRNA].keys(), 'clusters':miRNAs[miRNA], 'GO_terms':[]}
                        cc_go_miRNAs[miRNA]['GO_terms'].append(go1)


#################################################
## Filter network into cross-cancer subnetwork ##
#################################################

## Build the network: cancer -> miRNA -> cluster -> GO term -> hallmark of cancer
# Initialize directed graph
crossCancer = nx.DiGraph()

# Add each cancer -> miRNA -> cluster -> GO term -> hallmark of cancer relationship
for miRNA in cc_go_miRNAs:
    for cluster in cc_go_miRNAs[miRNA]['clusters']:
        # Add cancer -> miRNA
        crossCancer.add_edge(st8_hm.loc[cluster,'cancer'],miRDb.loc[miRNA,'Name'])
        crossCancer.nodes[st8_hm.loc[cluster,'cancer']]['group'] = 'cancer'
        crossCancer.nodes[miRDb.loc[miRNA,'Name']]['group'] = 'miRNA'
        # Add miRNA -> cluster
        crossCancer.add_edge(miRDb.loc[miRNA,'Name'],cluster)
        crossCancer.nodes[cluster]['group'] = 'cluster'
        # Add cluster -> GO term -> HM
        for go1 in cc_go_miRNAs[miRNA]['GO_terms']:
            crossCancer.add_edge(cluster,go1)
            crossCancer.nodes[go1]['group'] = 'go_term'
            for hm1 in cc_go_miRNAs[miRNA]['hallmarks']:
                crossCancer.add_edge(go1,hm1)
                crossCancer.nodes[hm1]['group'] = 'hallmark'

cc_nt = Network('800px','800px', directed=True)
cc_nt.show_buttons(filter_=['physics'])
cc_nt.from_nx(crossCancer)
cc_nt.show('crossCancer_wGOTerms.html')


## Build the network: cancer -> miRNA -> cluster -> hallmark of cancer
# Initialize directed graph
crossCancer = nx.DiGraph()

# Add each cancer -> miRNA -> cluster -> hallmark of cancer relationship
for miRNA in cc_go_miRNAs:
    for cluster in cc_go_miRNAs[miRNA]['clusters']:
        # Add cancer -> miRNA
        crossCancer.add_edge(st8_hm.loc[cluster,'cancer'],miRDb.loc[miRNA,'Name'])
        crossCancer.nodes[st8_hm.loc[cluster,'cancer']]['group'] = 'cancer'
        crossCancer.nodes[miRDb.loc[miRNA,'Name']]['group'] = 'miRNA'
        # Add miRNA -> cluster
        crossCancer.add_edge(miRDb.loc[miRNA,'Name'],cluster)
        crossCancer.nodes[cluster]['group'] = 'cluster'
        # Add cluster -> HM
        for hm1 in cc_go_miRNAs[miRNA]['hallmarks']:
            crossCancer.add_edge(cluster,hm1)
            crossCancer.nodes[hm1]['group'] = 'hallmark'

cc_nt = Network('800px','800px', directed=True)
cc_nt.show_buttons(filter_=['physics'])
cc_nt.from_nx(crossCancer)
cc_nt.show('crossCancer_no_GOTerms.html')


## Build the network: cancer -> miRNA -> cluster -> hallmark of cancer
# Initialize directed graph
crossCancer = nx.DiGraph()

# Add each cancer -> miRNA -> cluster -> hallmark of cancer relationship
for miRNA in cc_go_miRNAs:
    for cluster in cc_go_miRNAs[miRNA]['clusters']:
        # Add cancer -> miRNA
        crossCancer.add_edge(st8_hm.loc[cluster,'cancer'],miRDb.loc[miRNA,'Name'])
        crossCancer.nodes[st8_hm.loc[cluster,'cancer']]['group'] = 'cancer'
        crossCancer.nodes[miRDb.loc[miRNA,'Name']]['group'] = 'miRNA'
        # Add miRNA -> HM
        for hm1 in cc_go_miRNAs[miRNA]['hallmarks']:
            crossCancer.add_edge(miRDb.loc[miRNA,'Name'],hm1)
            crossCancer.nodes[hm1]['group'] = 'hallmark'

cc_nt = Network('800px','800px', directed=True)
cc_nt.show_buttons(filter_=['physics'])
cc_nt.from_nx(crossCancer)
cc_nt.show('crossCancer_no_Cluster_no_GOTerms.html')

