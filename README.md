# Proteomic_pipeline
R scripts for proteomic data analysis

The depostory contains three pipelines:

WGCNA:
Weighted correlation network analysis (WGCNA) is used for finding clusters (modules) of highly correlated genes/protiens. The central idea is that proteins that are correlated have some type of biological relatedness. Part of the pipeline, line 112, produces a table that is then used for GO-elite.

GO-elite:
The objective of GO-elite is to identify a  set of biological Ontology terms or pathways to describe a particular set of genes/proteins. It answeres the basic question who are these proteins and what do they do.

Random Forest:
The objective of this pipeline is to build a model that can predict diagnosis based of relative protein abundance. The pipeline also calculates feature importance. This information can be used to find biomarker candidates.

Together this pipelines can inform the researcher what group of proteins are related to disease. What pathways these proteins are in and what proteins are best used as biomarkers.

