# Proteomic_pipeline
R scripts for analysis proteimic data

The depostory contains two scripts one being WGCNA and the other being GO-elite

WGCNA
Weighted correlation network analysis (WGCNA) is used for finding clusters (modules) of highly correlated genes/protiens. The central idea is that proteins that are correlated have some type of biological relatedness. Part of the pipeline, line 112, produces a table that is then used for GO-elite.

GO-elite
The objective of GO-elite is to identify a  set of biological Ontology terms or pathways to describe a particular set of genes/proteins. It answeres the basic question who are these proteins and what do they do.

Together this pipeline can inform the researcher what group of proteins are related to disease. What pathways these proteins are.

