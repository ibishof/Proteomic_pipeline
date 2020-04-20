## One-Step GO-ELITE WITH USER PARAMETERS - by Eric Dammer, Divya Nandakumar
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
# 
# rm(list = ls())
dev.off()
 options(stringsAsFactors=FALSE)
# 
# ######################## EDIT THESE VARIABLES (USER PARAMETERS) #########################################################
fileName <- "test_inputs.csv"    #Input list, often comes from WGCNA pipeline
#             #INPUT is a CSV FILE - in the filePath folder.(The variable below)
#             #Can be formatted as Kme table from WGCNA pipeline, or
#             #can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
#             #in this case, the longest list is used as background for GO-Elite.
#             #  For simple columnwise list input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.
filePath <- "C:/Users/bishofij/Proteomics_Pipeline/go-elite_r_projects/"
#             #Folder that will contains the input file specified above, and where your output file will be created. 
#             #for system call from R to python to work, avoid folders with spaces; 
#             #but in Windows, spaces are handled by the script, making them passable.
# 
outFilename <- "test_inputs"
#             #Where you want the output to go. SUBFOLDER WITH THIS NAME WILL BE CREATED AUTOMATICALLY
# GOeliteFolder <- "C:/Users/bishofij/Proteomics_Pipeline/"
#             #has subfolder and python script for GO-Elite v1.2.5 (script should be edited per authors' instructions to Divya Nandakumar to correct bug for Fisher Exact Test
#             #Remove from end of file path GO-Elite_v.1.2.5-Win64/GO-Elite_v.1.2.5-Py
#             #for system call from R to python to work, avoid folders with spaces; 
#             #but in Windows, spaces are handled by the script, making them passable.
# 
 maxBarsPerOntology=5
#             #Ontologies per ontology type, used for generating the PDF report; does not limit GO-Elite output
# speciesCode="Hs"
#             #Hs for homo sapiens, Dm for fly; Mm, mouse; Rn, rat... (must have database downloaded via command line)
#             #if you use the GUI for GO-Elite, create a copy of the folders for command line, and delete databases downloaded via GUI.
# downloadDB=FALSE
#             #If TRUE, the database files for speciesCode species will be downloaded "(be patient)" from Ensembl (v62 preferred)
# pythonPath <- "C:/Python27/"
#             #python.exe for python v2.7 is here.
# 
panelDimensions=c(3,2)
#             #rows, columns for each page of PDF output
color=c("darkseagreen3","lightsteelblue1","lightpink4")
#             #colors respectively for ontologyTypes:
#             #"Biological Process","Molecular Function","Cellular Component"
#             #must be valid R colors
# modulesInMemory=FALSE
#             #uses cleanDat, net, and kMEdat from pipeline already in memory
# ANOVAgroups=FALSE
#             #if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
# ############ ADVANCED OPTION ####################################################################################
# customDBcmd=paste0("--customSet ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/Databases/EnsMart62Plus/C2/ --dataToAnalyze all ")
#             #set to "" if you have no custom database.
# customPanelDimensions=c(2,2)
# customReportMaxBars=20
# WGCNAinput=FALSE
# #Keep set to Flase
# ######################## REST OF THE CODE IS AUTOMATIC ##########################################################
# 
# 
# ## Clean out spaces and escaped backslashes from folder paths (folder names with spaces should not be used on non-windows systems with this script)
# filePath=paste0(paste( sapply(do.call(c,strsplit(filePath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
# pythonPath=paste0(paste( sapply(do.call(c,strsplit(pythonPath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
# GOeliteFolder=paste0(paste( sapply(do.call(c,strsplit(GOeliteFolder,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
# 
# 
# ## The input files for GO-Elite are text files with the gene list as the 1st column, a symbol identified (gene symbol, uniprot etc) as the 2nd column
# ## Different accepted inputs are given in the tutorial
# ## Commonly used symbols - Gene Symbol - Sy (example of input file below)
# ### GeneSymbol		SystemCode (Symbol format)
# ###	  GFAP		Sy
# ###	  APOE		Sy
# ## All input files are placed in one folder
# 
# ## The background file is prepared similarly and is placed in a separate folder
# ## The first part of the code runs GO-ELite either from R (using the system command) or can be run using the terminal (in mac)
# ## The second part call upon python and runs GO-elite. This requires GO-Elite to be installed and path to the GO-Elite installation site indicated following python.
# ## The 3rd part of the code plots the results from the GO-Elite results folder. This used are to create the plots.
# 
# ##-------------------------------##
# ## Preparing files for GO-Elite ##
# ## Takes in csv file breaks up the columns into seperate txt files
# 
# ##1.GO Elite input must be in a column-wise list format, and including longest such list as background.
# # We process the input file as simple lists by column in the CSV (largest list used as background)
# # Input and Denominator Files were prepared in the specified format (3 columns - genelist, format (Sy- for gene symbol) and kME)
# # Each module is a separate txt file. Denominator file contains background list whcih is all identified proteins in this case
# # All input files are placed in input or geneInfo folder and background file in denominator folder
# 
# #reread the file to a list of gene symbol (or UniqueID) lists
# dir.create(file.path(filePath, outFilename))
# 
#   #write lists to GOElite input files, and also the background file
#   for (i in names(DEXlistsForGO)) { 
#     dfGO<-data.frame(GeneSymbol=DEXlistsForGO[[i]],SystemCode=rep("Sy",length(DEXlistsForGO[[i]])))
#     write.table(unique(dfGO),file=paste(filePath,outFilename,"/",i,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
#   }
# 
# 
#   #write background
#   background <- unique(ANOVAout$Symbol)
#   background <- cbind(background,rep("Sy",length=length(background)))
#   colnames(background) <- c("GeneSymbol","SystemCode")
#   dir.create(file.path(paste0(filePath,outFilename),"background"))
#   write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
#   #nModules=length(names(DEXlistsForGO))
#   WGCNAinput=FALSE
# 
# 
#   modulesData <- as.list(read.csv(paste(filePath,fileName, sep=""),sep=",", stringsAsFactors=FALSE,header=T)) 
# 
#     nModules <- length(names(modulesData))
#   for (a in 1:nModules) {
#     modulesData[[a]] <- unique(modulesData[[a]][modulesData[[a]] != ""])
#     modulesData[[a]] <- modulesData[[a]][!is.na(modulesData[[a]])]
#     modulesData[[a]] <- do.call("rbind",strsplit(as.character(modulesData[[a]]), "[|]"))[,1]
#   }
#   ## Creating background file for GO Elite analysis
#   background <- modulesData[order(sapply(modulesData,length),decreasing=TRUE)][[1]]
#   background <- unique(background)
#   background <- cbind(background,rep("Sy",length=length(background)))
#   colnames(background) <- c("GeneSymbol","SystemCode")
#   dir.create(file.path(paste0(filePath,outFilename),"background"))
#   write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
# 
#   # Separate Symbol Lists into independent txt files for analysis by GO-Elite (CREATE INPUT FILES)
#   modulesData[[ names(modulesData[order(sapply(modulesData,length),decreasing=TRUE)])[1] ]] <- NULL
#   nModules = nModules -1 #no background
#   listNames <- uniquemodcolors <- names(modulesData)
#   for (i in listNames) {
#     listName <- i
#     listInfo <- cbind(modulesData[[listName]],rep("Sy",length=length(modulesData[[listName]])))
#     colnames(listInfo) <- c("GeneSymbol","SystemCode")
#     write.table(unique(listInfo),file=paste(filePath,outFilename,"/",listName,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
#   }
# 
# 
# 
# ##2. GO Elite Python Call
# ####----------------------- GO-ELite Analysis in Command Prompt/ Terminal or from R ------------------------------------#####
# 
# if (downloadDB) {
#   cat(paste0("Downloading Ensembl v62 database for species code ",speciesCode,". (Be patient.)\n\n"))
#   system( gsub("'",'\"',paste0(pythonPath,"python ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/GO_Elite.py --update Official --species ",speciesCode," --mod Ensembl --version EnsMart62Plus")) )
# }
# 
# commandLine=gsub("'",'"',paste0(pythonPath,"python ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/GO_Elite.py --species ",speciesCode," --mod Ensembl --permutations 'FisherExactTest' --method 'z-score' --zscore 1.96 --pval 0.05 --num 5 --input ",filePath,outFilename,"/ --denom ",filePath,outFilename,"/background/ ",customDBcmd,"--output ",filePath,outFilename,"/"))
# cat(paste0("NOW RUNNING THE FOLLOWING COMMAND:\n\n", commandLine,"\n\n(Estimated time for ", nModules, " lists to complete: ",round((30*nModules)/60,1)," minutes)\n Start time: ",Sys.time(),"\n"))
# system( gsub("'",'\"',paste0(pythonPath,"python ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/GO_Elite.py --species ",speciesCode," --mod Ensembl --permutations 'FisherExactTest' --method 'z-score' --zscore 1.96 --pval 0.05 --num 5 --input ",filePath,outFilename,"/ --denom ",filePath,outFilename,"/background/ ",customDBcmd,"--output ",filePath,outFilename,"/")) )



##3. Output Report of Z-Score Barplots, processing all GO-Elite output files
############################# ----------------------Plotting for modules ------------------------#######################
######## this script plots the top 3 ontologies for biological process, mol function and cell component for each module


##color scheme for ontology type key/legend (can be changed in user parameters, editing the "color" vector)
ontologyTypes=c("Biological Process","Molecular Function","Cellular Component")


modulesData <- as.list(read.csv(paste(filePath,fileName, sep=""),sep=",", stringsAsFactors=FALSE,header=T)) 

nModules <- length(names(modulesData))

# if(ANOVAgroups) {
#   xlabels <- names(DEXlistsForGO)
#   xlabels.frame <- data.frame(Colors=rep(NA,length(xlabels)),Labels=xlabels)
#   uniquemodcolors <- names(DEXlistsForGO) #not set above
# } else {
listNames <- uniquemodcolors <- names(modulesData)
xlabels <- uniquemodcolors #labels2colors(c(1:nModules))
xlabels1 <- paste("M",seq(1:nModules),sep="")
xlabels.frame <- as.data.frame(data.frame(Colors=xlabels,Labels=paste0(xlabels1," ",xlabels)))
#}

setwd(paste0(filePath,outFilename,"/"))
pdf(paste0("GO-Elite_",outFilename,".pdf"),height=10,width=8)
op <- par(mfrow=panelDimensions,oma=c(0,0,3,0))
frame()
legend(x="topleft",legend = ontologyTypes, fill=color, title=" ",cex=2,horiz=F,xpd=T)
legend(x="topleft",legend = c(" "," "," "), title="Ontology Types",cex=2.5,horiz=F,xpd=T, bty='n', title.adj=1.4)
#frame()
GOEliteOUTfileTrailer<- c("-GO_z-score_elite.txt")
summary <- list()
for(i in c(1:(length(uniquemodcolors)))){
  thismod=uniquemodcolors[i]	
  if (file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep="")) == F) { #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
    if(file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep="")) == T) {
      tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep=""),sep="\t")
    } else {
      next
    } 
  } else {
    tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep=""),sep="\t") #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
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

# 
# if(!customDBcmd=="") {
# library(stringr)
# pdf(paste0("GO-Elite_",outFilename,"-CUSTOM_db.pdf"),height=10,width=8)
# op <- par(mfrow=customPanelDimensions,oma=c(0,0,3,0))
# frame()
# legend(x="topleft",legend = ontologyTypes, fill=color, title=" ",cex=2,horiz=F,xpd=T)
# legend(x="topleft",legend = c(" "," "," "), title="Ontology Types",cex=2.5,horiz=F,xpd=T, bty='n', title.adj=1.4)
# 
# GOEliteOUTfileTrailer<-if(ANOVAgroups | !WGCNAinput) { c("-UserSuppliedAssociations_z-score_elite.txt"); } else { c("_Module-UserSuppliedAssociations_z-score_elite.txt"); }
# summary <- list()
# for(i in c(1:(length(uniquemodcolors)))){
# 	thismod=uniquemodcolors[i]	
# 	if (file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep="")) == F) { #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
# 		if(file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep="")) == T) {
# 			tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep=""),sep="\t")
# 		} else {
# 		next
# 		} 
# 	} else {
# 	tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep=""),sep="\t")
# 	}
# 	if (length(tmp[,2]) == 0) next
# 	tmp = tmp[,c(1,1,7,8,12)] ## Select GO-terms,GO-Type, Z-score,pValues and gene Lists
# 	tmp1 = tmp[order(tmp$Z.Score,decreasing=T),]
# 	tmp2 = tmp1
# 	tmp3 = tmp2
# 	tmp3 <- na.omit(tmp3)
# 	tmp3 <- tmp3[order(tmp3$Z.Score,decreasing=T),][c(1:customReportMaxBars),] #added this row, if you want to mix ontology types and sort by Z.Score only
# 	tmp3 <- tmp3[rev(rownames(tmp3)),]
# 
# 	summary[[i]] <- tmp3
# 	
# 	### To color bars by mol function, cell component or biological process
# 	for (j in 1:nrow(tmp3)){
# 		tmp3$color[j] <- color[1]
# 		if(WGCNAinput) { tmp3$color[j] <- uniquemodcolors[i] } #module color for all bars, instead of different colors by ontology type
# 	} 
# 
# #	if (tmp3$Z.Score == F) next
# 	if (is.na(max(tmp3$Z.Score))) tmp3<-na.omit(tmp3)
# 
# 	par(mar=c(4,15,4,3))
# 	xlim <- c(0,1.1*max(tmp3$Z.Score))	
# 	moduleTitle <- xlabels.frame[i,"Labels"]
# 	xh <- barplot(tmp3$Z.Score,horiz = TRUE,width =0.85,las=1,main=moduleTitle, xlim=xlim,col=tmp3$color,cex.axis=0.7,xlab="Z Score",cex.lab=0.9,cex.main=0.95,ylim=c(0,nrow(tmp3)+0.8))
# 	abline(v=1.96,col="red", cex.axis = 0.5)
# 	axis(2, at=xh, labels = str_to_title(gsub("_"," ",tmp3$Gene.Set.Name)), tick=FALSE, las =2, line =-0.5, cex.axis = 0.7) 
# }
# 
# par(op) # Leaves the last plot
# dev.off()
# } #end if(!customDBcmd=="")