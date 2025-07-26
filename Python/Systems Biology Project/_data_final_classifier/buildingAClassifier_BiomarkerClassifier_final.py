# -*- coding: utf-8 -*-
"""
@author: aditya
"""

#########################
## May need to install ##
#########################

#pip install GEOparse
#pip install fastcluster

####################
## Load libraries ##
####################

import GEOparse as gp
import pandas as pd
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, auc, roc_curve, RocCurveDisplay, confusion_matrix
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import ComplementNB
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

os.getcwd()
os.chdir("C:/Users/adi03/Documents/School/BME 524")

## Load up real-world data into pandas
#  - Training = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174302
#  - Validation = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142987
#gseTrain = gp.get_GEO(filepath="data_final_enrichment/GSE174302_family.soft/GSE174302_family.soft")
#gseValidation = gp.get_GEO(filepath="data_final_enrichment/GSE142987_family.soft/GSE174302_family.soft")

## Load up real-world data into pandas
#  - Training = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183795
##  - Validation = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28735
#  - Validation = https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452 ?maybe need this one for validation too
gseTrain = gp.get_GEO(filepath="data_final_enrichment/GSE183795_family.soft.gz")
gseValidation = gp.get_GEO(filepath="data_final_enrichment/GSE62452_family.soft.gz")

#pd.set_option('display.max_columns', None)
#pd.set_option('display.max_rows', None)

# GSE183795 (n = 244)
#gseTrain.phenotype_data.head()
phenosTrain = gseTrain.phenotype_data[['title','geo_accession','characteristics_ch1.0.tissue']]
phenosTrain.columns = ['title','geo_accession','tissue_type']
print(phenosTrain['tissue_type'].value_counts())

# GSE62452 (n = 130)
#gseValidation.phenotype_data.head()

phenosValidation = gseValidation.phenotype_data[['title','geo_accession','characteristics_ch1.0.tissue']]
phenosValidation.columns = ['title','geo_accession','tissue_type']
print(phenosValidation['tissue_type'].value_counts())

# Extract gene expression data
gexpTrain = pd.concat([gseTrain.gsms[gsm].table['VALUE'] for gsm in gseTrain.gsms.keys()],axis=1)
gexpTrain.index = list(gseTrain.gsms[list(gseTrain.gsms.keys())[0]].table['ID_REF'])
gexpTrain.columns = gseTrain.gsms.keys()
gexpValidation = pd.concat([gseValidation.gsms[gsm].table['VALUE'] for gsm in gseValidation.gsms.keys()],axis=1)
gexpValidation.index = list(gseValidation.gsms[list(gseValidation.gsms.keys())[0]].table['ID_REF'])
gexpValidation.columns = gseValidation.gsms.keys()

## Remove active TB data
gexpTrain = gexpTrain[phenosTrain.loc[phenosTrain['tissue_type']!='Normal pancreas'].index]
phenosTrain = phenosTrain.loc[phenosTrain['tissue_type']!='Normal pancreas']
#gexpValidation = gexpValidation[phenosValidation.loc[phenosValidation['tissue_type']!='Normal pancreas'].index]
#phenosValidation = phenosValidation.loc[phenosValidation['tissue_type']!='Normal pancreas']
print(phenosTrain['tissue_type'].value_counts())
print(phenosValidation['tissue_type'].value_counts())

## Feature selection is a very important step for classification.
# Too many noisy features can make detecting any signal difficult for most classification algroithms.
# Our feature selection will be selecting top most variant genes
# Goal is to remove invariant genes, which provide very little signal
# row are genes
top1000 = gexpTrain.var(axis=1).sort_values(ascending=False).index[range(1000)]
##########################
## Quatity control (QC) ##
## for data from GEO    ## Gene expression omnibus
##########################

## Quality control plots to ensure there are no outliers, and that sample labels are meaningful
# In PCA plot we want to see that similar disease states cluster, and that there aren't any outliers or batch effects.
# Clustermap is used in a similar way, but because their methods are different they are both useful.
with PdfPages('qc_GSE183795_Train150.pdf') as pdf:
    # QC using PCA anlaysis
    pca = PCA(n_components=2)
    gexp_pca = pca.fit(gexpTrain.loc[top1000].T).transform(gexpTrain.loc[top1000].T)
    plt.figure()
    colors = {'Tumor':'#fee5d9','adjacent non-tumor':'#0000ff'}
    for i in colors:
        plt.scatter(gexp_pca[phenosTrain['tissue_type']==i, 0], gexp_pca[phenosTrain['tissue_type']==i, 1],  color=colors[i], alpha=.8, linewidths=0.75, label=i, edgecolors='black')
    plt.legend(['Tumor','Adjacent Non-tumor'])
    plt.xlabel('PCA1')
    plt.ylabel('PCA2')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=False, scatterpoints=1)
    plt.title('PCA: Patients with Pancreatic Tumor vs\nAdjacent Non Tumor of pancreatic ductal adenocarcinoma')
    plt.tight_layout()
    pdf.savefig()
    plt.show()
    plt.close()
    # QC using heatmap to correlate samples based on top 2000 most variant genes
    # Note: in heatmap that not all axes labels are shown. But color bars are.
    c1 = gexpTrain.loc[top1000].corr()
    c1.columns = phenosTrain['tissue_type']
    c1.index = phenosTrain['tissue_type']
    colRow_colors = [colors[i] for i in phenosTrain['tissue_type']]
    sns.clustermap(c1, cmap='Blues', col_colors=colRow_colors, vmin=0, vmax=1, yticklabels=True,xticklabels=True)
    pdf.savefig()
    plt.show()
    plt.close() 


with PdfPages('qc_GSE62452_Validation150.pdf') as pdf:
    # QC using PCA anlaysis
    pca = PCA(n_components=2)
    gexp_pca = pca.fit(gexpValidation.loc[top1000].T).transform(gexpValidation.loc[top1000].T)
    plt.figure()
    
    colors = {'Pancreatic tumor':'#fee5d9','adjacent pancreatic non-tumor':'#0000ff'}
    for i in colors:
        plt.scatter(gexp_pca[phenosValidation['tissue_type']==i, 0], gexp_pca[phenosValidation['tissue_type']==i, 1], color=colors[i], alpha=.8, linewidths=0.75, label=i, edgecolors='black')
    plt.legend(['Tumor','Adjacent Non-tumor'])
    plt.xlabel('PCA1')
    plt.ylabel('PCA2')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), shadow=False, scatterpoints=1)
    plt.title('PCA: Patients with Pancreatic Tumor vs\nAdjacent Non Tumor of pancreatic ductal adenocarcinoma')
    plt.tight_layout()
    pdf.savefig()
    plt.show()
    plt.close()
    # QC using heatmap to correlate samples based on top 2000 most variant genes
    # Note: in heatmap that not all axes labels are shown. But color bars are.
    c1 = gexpValidation.loc[top1000].corr()
    c1.columns = phenosValidation['tissue_type']
    c1.index = phenosValidation['tissue_type']
    colRow_colors = [colors[i] for i in phenosValidation['tissue_type']]
    sns.clustermap(c1, cmap='Blues', col_colors=colRow_colors, vmin=0, vmax=1, yticklabels=True,xticklabels=True)
    pdf.savefig()
    plt.show()
    plt.close()
    

##########################################
## Building a classifier for Latent TB:        ##
## Classification with cross-validation ##
##########################################

## test dataset

gexpTrain110 = gexpTrain.loc[top1000].T
gexpTrain110_scld = StandardScaler().fit_transform(gexpTrain110)
gexpTrain110_scld_mm = MinMaxScaler().fit_transform(gexpTrain110_scld)

## Combine the disease state labels
disStatesTrain = phenosTrain['tissue_type']
print(disStatesTrain.value_counts())

## Generate labels that are in the set ['TB', 'other']
# This will bulid a classifier that can discriminate TB from the other normal and other lung diseases
disStatesTrain = pd.Series(['Tumor' if i=='Tumor' else 'other' for i in phenosTrain['tissue_type']])
print(disStatesTrain.value_counts())

## Set up classification comparison for 100 fold cross-validation
# Can change kFolds to change number of iterations
# Store out results in a list nested inside a dictionary
kFolds = 100
performanceResults = {'KNN':[],'SVM':[],'LDA':[],'RF':[],'CNB':[]}
for i in range(kFolds):
    ## Split up data into train and test sets using random integer seed
    gexpData_train, gexpData_test, disState_train, disState_test = train_test_split(gexpTrain110_scld, disStatesTrain, test_size=0.4, random_state=i)
    gexpData_train_mm, gexpData_test_mm, disState_train_mm, disState_test_mm = train_test_split(gexpTrain110_scld_mm, disStatesTrain, test_size=0.4, random_state=i)
    #
    #
    ## Train and test the classifiers
    # k-nearest neighbors (KNN)
    knn_clf = KNeighborsClassifier(8)
    knn_clf.fit(gexpData_train, disState_train)
    knn_pred = knn_clf.predict(gexpData_test)
    performanceResults['KNN'].append(classification_report(disState_test, knn_pred, output_dict=True, zero_division=0))
    # Support vector machines (SVM)
    svm_clf = SVC(kernel='linear', C = 0.025)
    svm_clf.fit(gexpData_train, disState_train)
    svm_pred = svm_clf.predict(gexpData_test)
    performanceResults['SVM'].append(classification_report(disState_test, svm_pred, output_dict=True, zero_division=0))
    # Linear discriminant analysis (LDA)
    lda_clf = LinearDiscriminantAnalysis()
    lda_clf.fit(gexpData_train, disState_train)
    lda_pred = lda_clf.predict(gexpData_test)
    performanceResults['LDA'].append(classification_report(disState_test, lda_pred, output_dict=True, zero_division=0))
    # Random Forest ensemble (RF)
    rf_clf = RandomForestClassifier()
    rf_clf.fit(gexpData_train, disState_train)
    rf_pred = rf_clf.predict(gexpData_test)
    performanceResults['RF'].append(classification_report(disState_test, rf_pred, output_dict=True, zero_division=0))
    # Complement Naive Bayes (CNB)
    cnb_clf = ComplementNB()
    cnb_clf.fit(gexpData_train_mm, disState_train)
    cnb_pred = cnb_clf.predict(gexpData_test_mm)
    performanceResults['CNB'].append(classification_report(disState_test, cnb_pred, output_dict=True, zero_division=0))


## Build figure to describe classifiers: using 'f1-score'
# f1 score = 2/(1/recall+1/precision) = tp/(tp + 0.5 * (fp + fn))
# tp = true positive, fp = false positive, tn = true negative, fn = false negative
f1_scores = {'KNN':{},'SVM':{},'LDA':{},'RF':{},'CNB':{}}
plotMe = pd.DataFrame(columns=['f1_score','method','class'])
for disState in list(set(disStatesTrain))+['weighted avg']: # ignoring macro average
    for clf in performanceResults.keys():
        f1_scores[clf][disState] = sum([performanceResults[clf][i][disState]['f1-score'] for i in range(kFolds)])/kFolds
        if not disState=='weighted avg':
            plotMe = pd.concat([plotMe,pd.DataFrame({'f1_score':[performanceResults[clf][i][disState]['f1-score'] for i in range(kFolds)], 'method':[clf for i in range(kFolds)], 'class':[disState for i in range(kFolds)]})])

# Print out the matrix of f1-scores for classifiers
print(pd.DataFrame(f1_scores))

# Make a plot of the f1-scores to compare across all 100 cross-validated f1-scores
sns.boxplot(x='class', y='f1_score', hue='method', data=plotMe)
plt.show()

## Calculate the positive predictive value (PPV) and negative predictive value (NPV)
# PPV = precision(TB) = tp / (tp + fp)
# NPV = prcision(other) = tn / (tn + fn)
# tp = true positive, fp = false positive, tn = true negative, fn = false negative
ppvNpv2 = {'KNN':{},'SVM':{},'LDA':{},'RF':{},'CNB':{}}
for disState in list(set(disStatesTrain))+['weighted avg']: # ignoring macro average
    for clf in performanceResults.keys():
        ppvNpv2[clf][disState] = sum([performanceResults[clf][i][disState]['precision'] for i in range(kFolds)])/kFolds

# Print out the matrix of postive and negative predictive power for classifiers pos for TB and neg for other
print(pd.DataFrame(ppvNpv2))

####################################################################
## Choose RF because:                                            ##
##  1. Best median f1-scores                                      ##
##  2. Best weighted avg precision(LTBI): RF PPV & precision(other):NPV ##
####################################################################

## Build final RF classifier by training on all combined data
# Train on full train/test dataset
rf_clf = RandomForestClassifier()
rf_clf.fit(gexpTrain110_scld, list(disStatesTrain))

## Plot most important genes in SVM classifier
# First make a DataFrame with all the information:  RF coefficient, probe ID and index by gene symbol (for human readability)
# GPL is the platform on the GEO website

importances = rf_clf.feature_importances_
std = np.std([tree.feature_importances_ for tree in rf_clf.estimators_], axis=0)
forest_importances = pd.Series(importances, index=gseTrain.gpls['GPL6244'].table['gene_assignment'].loc[gseTrain.gpls['GPL6244'].table['ID'].isin(gexpTrain110.columns)])
forest_importances_df = pd.DataFrame({'coef':forest_importances})
forest_importances_df['ID'] = forest_importances_df.index

forest_importances_df = forest_importances_df.loc[forest_importances_df['ID']!='---']
forest_importances_df['ID'] = forest_importances_df['ID'].str.split(' /// ')
forest_importances_df = forest_importances_df.explode('ID')
forest_importances_df['ID'] = forest_importances_df['ID'].str.split(' // ')

for i in range(0,len(forest_importances_df['ID'])):
    forest_importances_df['ID'][i] = forest_importances_df['ID'][i][1]
forest_importances_df = forest_importances_df.drop_duplicates()
forest_importances_df.index = forest_importances_df['ID']

# Make horizontal bar plot of top 25 genes
forest_importances_df = forest_importances_df.sort_values('coef', ascending=False)
forest_importances_df.iloc[list(reversed(list(range(0,25))))]['coef'].plot.barh()
plt.ylabel('RF Coefficient')
plt.tight_layout()
plt.show()
## Literature search of top genes:
# CD14 - found to be associated with TB/latent TB  (https://pmc.ncbi.nlm.nih.gov/articles/PMC8533229/)



##############
## Validate ##
##############

# Set up data for validation dataset
gexpValidation110 = gexpValidation.loc[top1000].T
gexpValidation110_scld = StandardScaler().fit_transform(gexpValidation110)
gexpValidation110_scld_mm = MinMaxScaler().fit_transform(gexpValidation110_scld)

# Set up disease states for validation dataset
disStatesValidation = pd.Series(['Tumor' if i=='Pancreatic tumor' else 'other' for i in phenosValidation['tissue_type']])
print(disStatesValidation.value_counts())

# Get results RF
rf_pred = rf_clf.predict(gexpValidation110_scld)
rf_pred_prob = rf_clf.predict_proba(gexpValidation110_scld)
print(pd.DataFrame(classification_report(disStatesValidation, rf_pred, output_dict=True, zero_division=0)))
print(confusion_matrix([int(i=='other') for i in disStatesValidation], [int(i=='other') for i in rf_pred]))

# Plot ROC curve for validation of RF
fig, ax = plt.subplots(figsize=(6, 6))
RocCurveDisplay.from_predictions([int(i=='other') for i in disStatesValidation], rf_pred_prob[:,1], name='RF', ax=ax)
plt.show()