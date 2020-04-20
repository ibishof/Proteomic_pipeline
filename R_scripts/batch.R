
BiocManager::install("limma")
library("limma")


BiocManager::install("sva")
library("sva")
library("limma")

#Installing and loading packages and data

BiocManager::install("bladderbatch")
library(bladderbatch)
data(bladderdata)

head.matrix(bladderEset)
#pdata(data, type, comment, metadata)
pheno = pData(bladderEset)

edata = exprs(bladderEset)

#Applying the ComBat function to adjust for known batches

#Pulls out batch data from pheno matrix
batch = pheno$batch

#Just the intercept value from pheno
modcombat = model.matrix(~1, data=pheno)

#The actual combat funcion. This outputs the batch corrected data set "combat_edata)
combat_edata = ComBat(dat=edata, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE)
