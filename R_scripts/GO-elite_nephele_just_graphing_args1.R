
# ##############################################################################
# ## One-Step GO-ELITE WITH USER PARAMETERS - by Eric Dammer, Divya Nandakumar
# ##  - performs GO-Elite v1.2.5 with Fisher Exact Test for enrichment p<0.05
# ##    and 5 minimum genes per ontology
# ##
# ## Nicholas Seyfried Lab Bioinformatics - for the lab - 02/11/2019 update
# ##  - Change parameters between the lines below
# ##############################################################################
# ### GO-Elite is a python package to perform ontology analysis on gene sets. 
# ### The python script for GO-Elite can be downloaded from http://www.genmapp.org/go_elite/
# ### Alternatively, there is also a GUI which can be downlaoded from the same website.
# ###
# ### Custom databases can be downloaded from http://software.broadinstitute.org/gsea/msigdb/
# ### including TRANSFAC and curated pathways collected from various sources, e.g. the C2 DB.
# ### GO-Elite requires python 2.7 be installed, and FET with command-line requires
# ### that a bugfix be applied; copy GO_Elite.py included to GO-Elite program subfolder:
# ### GOeliteFolder/GO-Elite_v.1.2.5-Py/
# ##############################################################################
# 
# 
########################## Prepareing Enviromeent for Graphing #################################################
 rm(list = ls()) # Removed objects in the R enviroment. Not sure if needed on AWI
 dev.off() # Turns off graphic devices to allow this intance of the code to use graphing later. Not sure if needed on AWI
 dev.off()
 options(stringsAsFactors=FALSE)
 options(error = browser())
 traceback()
# 
# ######################## EDIT THESE VARIABLES (USER PARAMETERS) #########################################################

args <- c("test_inputs.csv", file.path("C:", "Users", "bishofij", "Proteomics_Pipeline", "go-elite_r_projects"),"test_inputs", file.path("GO-Elite_results", "CompleteResults","ORA_pruned"))
 
go_elite_output <- file.path("GO-Elite_results", "CompleteResults","ORA_pruned") 
              #This is the path to the GO-elite statistical results
#########################################################################################
# 
maxBarsPerOntology=5
#             #Ontologies per ontology type, used for generating the PDF report; does not limit GO-Elite output
panelDimensions=c(3,2)
#             #rows, columns for each page of PDF output
color=c("darkseagreen3","lightsteelblue1","lightpink4")
#             #colors respectively for ontologyTypes:
#             #"Biological Process","Molecular Function","Cellular Component"
#             #must be valid R colors
ontologyTypes=c("Biological Process","Molecular Function","Cellular Component")

##3. Output Report of Z-Score Barplots, processing all GO-Elite output files
############################# ----------------------Plotting for modules ------------------------#######################
######## this script plots the top 3 ontologies for biological process, mol function and cell component for each module


##This set up the columns that will be used for the graphing

modulesData <- as.list(read.csv(file.path(args[2],args[1]), stringsAsFactors=FALSE,header=T))
listNames <- uniquemodcolors <- names(modulesData)
nModules <- length(names(modulesData))

xlabels <- uniquemodcolors #labels2colors(c(1:nModules))
xlabels1 <- paste("M",seq(1:nModules),sep="")
xlabels.frame <- as.data.frame(data.frame(Colors=xlabels,Labels=paste0(xlabels1," ",xlabels)))
#}

setwd(paste(file.path(args[2],args[3])))
pdf(paste0("GO-Elite_",args[3],".pdf"),height=10,width=8)
op <- par(mfrow=panelDimensions,oma=c(0,0,3,0))
frame()
legend(x="topleft",legend = ontologyTypes, fill=color, title=" ",cex=2,horiz=F,xpd=T)
legend(x="topleft",legend = c(" "," "," "), title="Ontology Types",cex=2.5,horiz=F,xpd=T, bty='n', title.adj=1.4)
#frame()
GOEliteOUTfileTrailer<- c("-GO_z-score_elite.txt")
summary <- list()

############################################# The Visualization #####################################################
for(i in c(1:(length(uniquemodcolors)))){
  thismod=uniquemodcolors[i]	
  if (file.exists(file.path(args[2],args[3],go_elite_output,paste(thismod,GOEliteOUTfileTrailer, sep=""))) == F) { #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
    if(file.exists(file.path(args[2],args[3],go_elite_output, paste(thismod,"_2",GOEliteOUTfileTrailer, sep=""))) == T) {
      tmp=read.csv(file=file.path(args[2],args[3],go_elite_output,paste(thismod,"_2",GOEliteOUTfileTrailer, sep="")),sep="\t")
    } else {
      next
    } 
  } else {
    tmp=read.csv(file=file.path(args[2],args[3],go_elite_output,paste(thismod,GOEliteOUTfileTrailer, sep="")),sep="\t") #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
  }
  if (length(tmp[,2]) == 0) next
  tmp = tmp[,c(2,3,9,10,11)] ## Select GO-terms,GO-Type, Z-score,pValues and gene Lists
  tmp1 = tmp[order(tmp$Z.Score,decreasing=T),]
  tmp2 = tmp1[order(tmp1$Ontology.Type,decreasing=T),] #was tmp2
  tmp3 = tmp2[tmp2$Ontology.Type == "biological_process",][c(1:maxBarsPerOntology),]
  tmp3 = rbind(tmp3,tmp2[tmp2$Ontology.Type == "molecular_function",][c(1:maxBarsPerOntology),] )
  tmp3 = rbind(tmp3,tmp2[tmp2$Ontology.Type == "cellular_component",][c(1:maxBarsPerOntology),] )
  tmp3 <- na.omit(tmp3)
  #	tmp3 <- tmp3[order(tmp3$Z.Score,decreasing=T),] #added this row, if you want to mix ontology types and sort by Z.Score only
  tmp3 <- tmp3[rev(rownames(tmp3)),]
  
  summary[[i]] <- tmp3
  
  ### To color bars by mol function, cell component or biological process
  for (j in 1:nrow(tmp3)){
    if (tmp3$Ontology.Type[j] == "molecular_function"){
      tmp3$color[j] <- color[2]
    } else if (tmp3$Ontology.Type[j] == "cellular_component"){
      tmp3$color[j] <- color[3]
    } else if (tmp3$Ontology.Type[j] == "biological_process"){
      tmp3$color[j] <- color[1]
    }
    # tmp3$color[j] <- uniquemodcolors[i] #module color for all bars, instead of different colors by ontology type
  } 
  
  if (tmp3$Z.Score == F) next
  par(mar=c(4,15,4,3))
  xlim <- c(0,1.1*max(tmp3$Z.Score))	
  moduleTitle <- xlabels.frame[i,"Labels"]
  xh <- barplot(tmp3$Z.Score,horiz = TRUE,width =0.85,las=1,main=moduleTitle, xlim=xlim,col=tmp3$color,cex.axis=0.7,xlab="Z Score",cex.lab=0.9,cex.main=0.95,ylim=c(0,nrow(tmp3)+0.8))
  abline(v=1.96,col="red", cex.axis = 0.5)
  axis(2, at=xh, labels = tmp3$Ontology.Name, tick=FALSE, las =2, line =-0.5, cex.axis = 0.7) 
}

par(op) # Leaves the last plot
dev.off()

