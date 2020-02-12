###MetabolomicsParser
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""This module contains methods for reading the HMDB and storing relationships"""

import sys, string
import os.path
import unique
import export
import time
import update; reload(update)
import OBO_import
import gene_associations
import traceback

############# Common file handling routines ############# 
def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv": dir_list2.append(entry)
    return dir_list2

def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

def lowerSymbolDB(source_to_gene):
    source_to_gene2={}
    for symbol in source_to_gene:
        source_to_gene2[string.lower(symbol)]=source_to_gene[symbol]
    return source_to_gene2
    
def verifyFile(filename):
    fn=filepath(filename); file_found = 'yes'
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception: file_found = 'no'
    return file_found

def importSpeciesData():
    if program_type == 'GO-Elite': filename = 'Config/species_all.txt' ### species.txt can be cleared during updating
    else: filename = 'Config/goelite_species.txt'
    x=0
    fn=filepath(filename);global species_list; species_list=[]; global species_codes; species_codes={}
    global species_names; species_names={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t'); abrev=t[0]; species=t[1]
        if x==0: x=1
        else:
            species_list.append(species)
            species_codes[species] = abrev
            species_names[abrev] = species

def getSourceData():
    filename = 'Config/source_data.txt'; x=0
    fn=filepath(filename)
    global source_types; source_types={}
    global system_codes; system_codes={}
    global mod_types; mod_types=[]
    for line in open(fn,'rU').readlines():
        data = cleanUpLine(line)
        t = string.split(data,'\t'); source=t[0]
        try: system_code=t[1]
        except IndexError: system_code = 'NuLL'
        if x==0: x=1
        else:
            if len(t)>2: ### Therefore, this ID system is a potential MOD
                if t[2] == 'MOD': mod_types.append(source)
            source_types[source]=system_code
            system_codes[system_code] = source ###Used when users include system code data in their input file
            
############# File download/extraction #############           
def downloadPAZARAssocations():
    base_url = 'http://www.pazar.info/tftargets/'
    filenames = getPAZARFileNames()
    print 'Downloading Transcription Factor to Target associations'
    source = 'raw'
    r = 4; k = -1
    for resource in filenames:
        filename = filenames[resource]
        url = base_url+filename
        start_time = time.time()
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/PAZAR/','')
        end_time = time.time()
        if (end_time-start_time)>3: ### Hence the internet connection is very slow (will take forever to get everything)
            downloadPreCompiledPAZAR() ### Just get the compiled symbol data instead
            print '...access to source PAZAR files too slow, getting pre-compiled from genmapp.org'
            source = 'precompiled'
            break
        k+=1
        if r==k:
            k=0
            print '*',
    print ''
    return source

def downloadPreCompiledPAZAR():
    """ Downloads the already merged symbol to TF file from PAZAR files """
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/tf-target.txt'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/PAZAR/symbol/','')
     
def downloadAmadeusPredictions():
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/symbol-Metazoan-Amadeus.txt'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Amadeus/','')
    
def downloadBioMarkers():
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/Hs_exon_tissue-specific_protein_coding.zip'
    print 'Downloading BioMarker associations'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/BioMarkers/','')
    
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/Mm_gene_tissue-specific_protein_coding.zip'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/BioMarkers/','')
    
def downloadKEGGPathways(species):
    print "Integrating KEGG associations for "+species
    url = 'http://www.genmapp.org/go_elite/Databases/KEGG/'+species+'-KEGG_20110518.zip'
    ### This is a fixed date resource since KEGG licensed their material after this date
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/KEGG/','')

def downloadDomainAssociations(selected_species):
    paths=[]
    if selected_species != None: ### Restrict to selected species only
        current_species_dirs=selected_species
    else:
        current_species_dirs = unique.read_directory('/'+database_dir)
    for species in current_species_dirs:
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/Domains/'+species+'_Ensembl-Domain.gz'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Domains/','txt')
        if 'Internet' not in status:
            paths.append((species,fln))
    return paths

def downloadPhenotypeOntologyOBO():
    print 'Downloading Phenotype Ontology structure and associations'
    url = 'ftp://ftp.informatics.jax.org/pub/reports/MPheno_OBO.ontology'
    fln,status = update.downloadSuppressPrintOuts(url,'OBO/','')

def downloadPhenotypeOntologyGeneAssociations():
    url = 'ftp://ftp.informatics.jax.org/pub/reports/HMD_HumanPhenotype.rpt'
    #url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/HMD_HumanPhenotype.rpt'
    ### Mouse and human gene symbols and gene IDs (use the gene symbols)
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Pheno/','')

def downloadPathwayCommons():
    print 'Downloading PathwayCommons associations'
    url = 'http://www.pathwaycommons.org/pc-snapshot/current-release/gsea/by_species/homo-sapiens-9606-gene-symbol.gmt.zip'
    fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/PathwayCommons/','')
    
def downloadDiseaseOntologyOBO():
    print 'Downloading Disease Ontology structure and associations'
    
    """ Unfortunately, we have to download versions that are not as frequently updated, since RGDs server
        reliability is poor """
    #url = 'ftp://rgd.mcw.edu/pub/data_release/ontology_obo_files/disease/CTD.obo'
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/CTD.obo'
    
    ### Includes congenital and environmental diseases - http://ctdbase.org/detail.go?type=disease&acc=MESH%3aD002318
    fln,status = update.downloadSuppressPrintOuts(url,'OBO/','')
    
def downloadDiseaseOntologyGeneAssociations(selected_species):
    if selected_species == None: sc = []
    else: sc = selected_species
    """ Unfortunately, we have to download versions that are not as frequently updated, since RGDs server
        reliability is poor """
        
    if 'Hs' in sc or len(sc)==0:
        #url = 'ftp://rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/homo_genes_do'
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/homo_genes_do'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Disease/','')

    if 'Mm' in sc or len(sc)==0:
        #url = 'ftp://rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/mus_genes_do'
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/mus_genes_do'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Disease/','')

    if 'Rn' in sc or len(sc)==0:
        #url = 'ftp://rgd.mcw.edu/pub/data_release/annotated_rgd_objects_by_ontology/rattus_genes_do'
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/rattus_genes_do'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/Disease/','')

def downloadMiRDatabases(species):
    url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/'+species+'_microRNA-Ensembl-GOElite_strict.txt'
    selected = ['Hs','Mm','Rn'] ### these are simply zipped where the others are not
    ### These files should be updated on a regular basis
    if species in selected:
        url = string.replace(url,'.txt','.zip')
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/microRNATargets/','')
    else:
        ### Where strict is too strict
        url = 'http://www.genmapp.org/go_elite/Databases/ExternalSystems/'+species+'_microRNA-Ensembl-GOElite_lax.txt'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/microRNATargets/','')
    fln = string.replace(fln,'.zip','.txt')
    return fln

def downloadGOSlimOBO():
    url = 'http://www.geneontology.org/GO_slims/goslim_pir.obo'
    #url = 'http://www.geneontology.org/GO_slims/goslim_generic.obo'  ### Missing 
    fln,status = update.downloadSuppressPrintOuts(url,'OBO/','')
    
def importUniProtAnnotations(species_db):
    base_url = 'http://www.altanalyze.org/archiveDBs/'
    uniprot_ensembl_db={}
    for species in species_db:
        url = base_url+species+'/custom_annotations.txt'
        fln,status = update.downloadSuppressPrintOuts(url,'BuildDBs/UniProt/'+species+'/','')
        for line in open(fln,'rU').xreadlines():
            data = cleanUpLine(line)
            ens_gene,compartment,function,symbols,full_name,uniprot_name,uniprot_ids,unigene = string.split(data,'\t')
            symbols = string.split(string.replace(symbols,'; Synonyms=',', '),', ')
            uniprot_ensembl_db[species,uniprot_name] = ens_gene
            species_extension = string.split(uniprot_name,'_')[-1]
            full_name = string.split(full_name,';')[0]
            if 'Transcription factor' in full_name:
                symbols.append(string.split(full_name,'Transcription factor ')[-1]) ### Add this additional synonym to symbols
            ### Extend this database out to account for weird names in PAZAR
            for symbol in symbols:
                new_name = string.upper(symbol)+'_'+species_extension
                if new_name not in uniprot_ensembl_db:
                    uniprot_ensembl_db[species,symbol+'_'+species_extension] = ens_gene
                uniprot_ensembl_db[species,string.upper(symbol)] = ens_gene
    return uniprot_ensembl_db

############# Import/processing/export #############    
def getPAZARFileNames():
    """ Filenames are manually and periodically downloaded from: http://www.pazar.info/cgi-bin/downloads_csv.pl"""
    fn = filepath('Config/PAZAR_list.txt')
    x=0
    filenames = {}
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            resource, filename = string.split(data,'\t')
            filenames[resource]=filename
    return filenames

class TFTargetInfo:
    def __init__(self,tf_name,ens_gene,project,pmid,analysis_method):
        self.tf_name=tf_name
        self.ens_gene=ens_gene
        self.project=project
        self.pmid=pmid
        self.analysis_method=analysis_method
    def TFName(self): return self.tf_name
    def Ensembl(self): return self.ens_gene
    def Project(self):
        if self.project[-1]=='_':
            return self.project[:-1]
        else:
            return self.project
    def PMID(self): return self.pmid
    def AnalysisMethod(self): return self.analysis_method
    def __repr__(self): return self.TFName()
    
def importPAZARAssociations():
    pazar_files = unique.read_directory('/BuildDBs/PAZAR')
    species_db={}
    tf_to_target={}
    for file in pazar_files:
        if '.csv' in file:
            name = string.join(string.split(file,'_')[1:-1],'_')
            fn = filepath('BuildDBs/PAZAR/'+file)
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                try:
                    ### Each line contains the following 11 tab-delim fields:
                    ### Fields are: <PAZAR TF ID>  <TF Name>  <PAZAR Gene ID>  <ensembl gene accession>  <chromosome>  <gene start coordinate>  <gene end coordinate>  <species>  <project name>  <PMID>  <analysis method> 
                    pazar_tf_id, tf_name, pazar_geneid, ens_gene, chr, gene_start,gene_end,species,project,pmid,analysis_method = string.split(data,'\t')
                    species,genus = string.split(species,' ')
                    species = species[0]+genus[0]
                    tft=TFTargetInfo(tf_name,ens_gene,project,pmid,analysis_method)
                    try: tf_to_target[species,tf_name].append(tft)
                    except Exception: tf_to_target[species,tf_name] = [tft]
                    species_db[species]=[]
                except Exception:
                    None ### Occurs due to file formatting issues (during an update?)

    determine_tf_geneids = 'no'
    if determine_tf_geneids == 'yes':
        """ The below code is probably most useful for creation of complex regulatory inference networks in Cytoscape """
        uniprot_ensembl_db = importUniProtAnnotations(species_db)
        missing=[]
        tf_to_target_ens={}
        for (species,tf_name) in tf_to_target:
            original_tf_name = tf_name
            try:
                ens_gene = uniprot_ensembl_db[species,tf_name]
                tf_to_target_ens[ens_gene]=tf_to_target[species,tf_name]
            except Exception:
                try:
                    tf_name = string.split(tf_name,'_')[0]
                    ens_gene = uniprot_ensembl_db[species,tf_name]
                    tf_to_target_ens[ens_gene]=tf_to_target[species,original_tf_name]
                except Exception:
                    try:
                        tf_names=[]
                        if '/' in tf_name:
                            tf_names = string.split(tf_name,'/')
                        elif ' ' in tf_name:
                            tf_names = string.split(tf_name,' ')
                        for tf_name in tf_names:
                            ens_gene = uniprot_ensembl_db[species,tf_name]
                            tf_to_target_ens[ens_gene]=tf_to_target[species,original_tf_name]          
                    except Exception: missing.append((tf_name,species))
        print 'Ensembl IDs found for UniProt Transcription factor names:',len(tf_to_target_ens),'and missing:', len(missing)
        #print missing[:20]
        
    ### Translate all species data to gene symbol to export for all species
    species_tf_targets={}
    for (species,tf_name) in tf_to_target:
        try:
            tf_db = species_tf_targets[species]
            tf_db[tf_name] = tf_to_target[species,tf_name]
        except Exception:
            tf_db = {}
            tf_db[tf_name] = tf_to_target[species,tf_name]
            species_tf_targets[species] = tf_db
        
    tf_dir = 'BuildDBs/PAZAR/symbol/tf-target.txt'
    tf_data = export.ExportFile(tf_dir)
    tf_to_symbol={}
    #print 'Exporting:',tf_dir
    #print len(species_tf_targets)
    for species in species_tf_targets:
        try: gene_to_source_id = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
        except Exception: gene_to_source_id={}
        tf_db = species_tf_targets[species]
        for tf_name in tf_db:
            for tft in tf_db[tf_name]:
                try:
                    for symbol in gene_to_source_id[tft.Ensembl()]:
                        symbol = string.lower(symbol)
                        tf_id = tf_name+'(Source:'+tft.Project()+'-PAZAR'+')'
                        tf_data.write(tf_id+'\t'+symbol+'\n')
                        try: tf_to_symbol[tf_id].append(symbol)
                        except Exception: tf_to_symbol[tf_id] = [symbol]
                except Exception: null=[]; 
    tf_data.close()
    tf_to_symbol = gene_associations.eliminate_redundant_dict_values(tf_to_symbol)
    return tf_to_symbol

def importPAZARcompiled():
    """ Skips over the above function when these tf-target file is downlaoded directly """
    tf_dir = 'BuildDBs/PAZAR/symbol/tf-target.txt'
    tf_to_symbol={}
    fn = filepath(tf_dir)
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        tf_id,symbol = string.split(data,'\t')
        try: tf_to_symbol[tf_id].append(symbol)
        except Exception: tf_to_symbol[tf_id] = [symbol]
    tf_to_symbol = gene_associations.eliminate_redundant_dict_values(tf_to_symbol)
    return tf_to_symbol

def importPhenotypeOntologyGeneAssociations():
    x=0
    pheno_symbol={}; phen=[]
    fn = filepath('BuildDBs/Pheno/HMD_HumanPhenotype.rpt')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            t = string.split(data,'\t')
            hs_symbol=t[0]; hs_entrez=t[1]; mm_symbol=t[2]; mgi=t[3]; pheno_ids=t[4]
            hs_symbol = string.lower(hs_symbol)
            mm_symbol = string.lower(mm_symbol)
            symbols = [mm_symbol,hs_symbol]
            pheno_ids = string.split(pheno_ids,' '); phen+=pheno_ids
            for pheno_id in pheno_ids:
                if len(pheno_id)>0:
                    for symbol in symbols:
                        try: pheno_symbol[pheno_id].append(symbol)
                        except Exception: pheno_symbol[pheno_id]=[symbol]
    phen = unique.unique(phen)
    pheno_symbol = gene_associations.eliminate_redundant_dict_values(pheno_symbol)
    return pheno_symbol

def importAmandeusPredictions(force):
    if force == 'yes':
        downloadAmadeusPredictions()
        
    x=0
    tf_symbol_db={}
    fn = filepath('BuildDBs/Amadeus/symbol-Metazoan-Amadeus.txt')
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            symbol,system,tf_name = string.split(data,'\t')
            symbol = string.lower(symbol)
            try: tf_symbol_db[tf_name].append(symbol)
            except Exception: tf_symbol_db[tf_name]=[symbol]
    tf_symbol_db = gene_associations.eliminate_redundant_dict_values(tf_symbol_db)
    return tf_symbol_db
    
def importDiseaseOntologyGeneAssocations():
    disease_ontology_files = unique.read_directory('/BuildDBs/Disease')
    symbol_to_DO={}
    for file in disease_ontology_files:
        if '_do' in file:
            fn = filepath('BuildDBs/Disease/'+file)
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if len(t)>1:
                    symbol=string.lower(t[2]); doid = t[4]
                    try: symbol_to_DO[doid].append(symbol)
                    except Exception: symbol_to_DO[doid]=[symbol]
    return symbol_to_DO
    
def exportSymbolRelationships(pathway_to_symbol,selected_species,pathway_type,type):    
    if selected_species != None: ### Restrict to selected species only
        current_species_dirs=selected_species
    else:
        current_species_dirs = unique.read_directory('/'+database_dir)
    
    for species in current_species_dirs:
        if '.' not in species:
            ens_dir = database_dir+'/'+species+'/gene-'+type+'/Ensembl-'+pathway_type+'.txt'
            ens_data = export.ExportFile(ens_dir)
            if 'mapp' in type: ens_data.write('GeneID\tSystem\tGeneSet\n')
            else: ens_data.write('GeneID\tGeneSet\n')
            try: ens_to_entrez = gene_associations.getGeneToUid(species,('hide','Ensembl-EntrezGene'))
            except Exception: ens_to_entrez ={}
            if len(ens_to_entrez)>0:
                entrez_dir = database_dir+'/'+species+'/gene-'+type+'/EntrezGene-'+pathway_type+'.txt'
                entrez_data = export.ExportFile(entrez_dir)
                if 'mapp' in type: entrez_data.write('GeneID\tSystem\tGeneSet\n')
                else: entrez_data.write('GeneID\tGeneSet\n')
            #print 'Exporting '+pathway_type+' databases for:',species
            try: gene_to_source_id = gene_associations.getGeneToUid(species,('hide','Ensembl-Symbol'))
            except Exception: gene_to_source_id={}
            source_to_gene = OBO_import.swapKeyValues(gene_to_source_id)
            source_to_gene = lowerSymbolDB(source_to_gene)
            for pathway in pathway_to_symbol:
                for symbol in pathway_to_symbol[pathway]:
                    try:
                        genes = source_to_gene[symbol]
                        for gene in genes:
                            if 'mapp' in type: ens_data.write(gene+'\tEn\t'+pathway+'\n')
                            else: ens_data.write(gene+'\t'+pathway+'\n')
                            if gene in ens_to_entrez:
                                for entrez in ens_to_entrez[gene]:
                                    if 'mapp' in type: entrez_data.write(entrez+'\tL\t'+pathway+'\n')
                                    else: entrez_data.write(entrez+'\t'+pathway+'\n')
                    except Exception: null=[]
            ens_data.close()
            try: entrez_data.close()
            except Exception: null=[]
            
def extractKEGGAssociations(species,mod,system_codes):
    import_dir = filepath('/BuildDBs/KEGG')
    g = gene_associations.GrabFiles(); g.setdirectory(import_dir)
    filedir = g.getMatchingFolders(species)
    gpml_data,pathway_db = gene_associations.parseGPML(filepath(filedir))
    gene_to_WP = gene_associations.unifyGeneSystems(gpml_data,species,mod)

    gene_associations.exportCustomPathwayMappings(gene_to_WP,mod,system_codes,filepath(database_dir+'/'+species+'/gene-mapp/'+mod+'-KEGG.txt'))

def extractGMTAssociations(species,mod,system_codes,data_type):
    if mod != 'HMDB':
        import_dir = filepath('/BuildDBs/'+data_type)
        gmt_data = gene_associations.parseGMT(import_dir)
        gene_to_custom = gene_associations.unifyGeneSystems(gmt_data,species,mod)
        gene_associations.exportCustomPathwayMappings(gene_to_custom,mod,system_codes,filepath(database_dir+'/'+species+'/gene-mapp/'+mod+'-'+data_type+'.txt'))
        
def transferGOSlimGeneAssociations(selected_species):
    if selected_species != None: ### Restrict to selected species only
        current_species_dirs=selected_species
    else:
        current_species_dirs = unique.read_directory('/'+database_dir)
    for species_code in current_species_dirs:
        try:
            ens_go_file_dir = filepath(database_dir+'/'+species_code+'/gene-go/Ensembl-GOSlim.txt')
            goslim_ens_file = filepath(database_dir+'/'+species_code+'/uid-gene/Ensembl-goslim_goa.txt')
            export.copyFile(goslim_ens_file,ens_go_file_dir)
            translateToEntrezGene(species_code,ens_go_file_dir)
        except Exception: null=[]

def translateToEntrezGene(species,filename):
    x=0; type = 'pathway'
    try: ens_to_entrez = gene_associations.getGeneToUid(species,('hide','Ensembl-EntrezGene'))
    except Exception: ens_to_entrez ={}

    if len(ens_to_entrez)>0:
        export_file = string.replace(filename,'Ensembl','EntrezGene')
        export_data = export.ExportFile(export_file)
        export_data.write('EntrezGene\tOntologyID\n')
        fn = filepath(filename)
        for line in open(fn,'rU').xreadlines():
            if x==0: x=1
            else:
                data = cleanUpLine(line)
                try:
                    ensembl,pathway = string.split(data,'\t')
                    type = 'ontology'
                except Exception:
                    ensembl,null,pathway = string.split(data,'\t')
                try:
                    entrezs = ens_to_entrez[ensembl]
                    for entrez in entrezs:
                        if type == 'ontology':
                            export_data.write(entrez+'\t'+pathway+'\n')
                        else:
                            export_data.write(entrez+'\tEn\t'+pathway+'\n')
                except Exception:
                    null=[]
        export_data.close()

def importMiRGeneAssociations(species_code,source_path):
    try:
        destination_path = filepath(database_dir+'/'+species_code+'/gene-mapp/Ensembl-microRNATargets.txt')
        export.copyFile(source_path,destination_path)
        translateToEntrezGene(species_code,destination_path)
    except Exception: null=[]            
        
def importBioMarkerGeneAssociations():
    try:
        biomarker_files = unique.read_directory('BuildDBs/BioMarkers/')
    except Exception:
        biomarker_files = unique.read_directory('/BuildDBs/BioMarkers/')
    x=0; marker_symbol_db={}
    for file in biomarker_files:
        if '.txt' in file:
            fn = filepath('BuildDBs/BioMarkers/'+file)
            for line in open(fn,'rU').xreadlines():
                data = cleanUpLine(line)
                t = string.split(data,'\t')
                if x==0:
                    x = 1; y=0
                    for i in t:
                        if 'marker-in' in i: mi = y
                        if 'Symbol' in i: sy = y
                        y+=1
                ensembl = t[0]; symbol = string.lower(t[sy]); marker = t[mi]
                markers = string.split(marker,'|')
                for marker in markers:
                    try: marker_symbol_db[marker].append(symbol)
                    except Exception: marker_symbol_db[marker]=[symbol]
    marker_symbol_db = gene_associations.eliminate_redundant_dict_values(marker_symbol_db)
    return marker_symbol_db

def importDomainGeneAssociations(species_code,source_path):
    try:
        destination_path = filepath(database_dir+'/'+species_code+'/gene-mapp/Ensembl-Domains.txt')
        export.copyFile(source_path,destination_path)
        translateToEntrezGene(species_code,destination_path)
    except Exception: null=[]            
          
############# Central buid functions #############

def importWikiPathways(selected_species,force):
    if selected_species == None:
        selected_species = unique.read_directory('/'+database_dir)
    importSpeciesData()
    getSourceData()
    all_species = 'no'
    if force == 'yes':
        try:
            gene_associations.convertAllGPML(selected_species,all_species) ### Downloads GPMLs and builds flat files
            status = 'built'
        except IOError:
            print 'Unable to connect to http://www.wikipathways.org'
            status = 'failed'
    status = 'built'
    if status == 'built':
        import BuildAffymetrixAssociations

        for species_code in selected_species:
            species_name = species_names[species_code]
            if status == 'built':          
                relationship_types = ['native','mapped']
                for relationship_type in relationship_types:
                    #print 'Processing',relationship_type,'relationships'
                    index=0
                    integrate_affy_associations = 'no'
                    incorporate_previous = 'yes'
                    process_affygo = 'no'
                    counts = BuildAffymetrixAssociations.importWikipathways(source_types,incorporate_previous,process_affygo,species_name,species_code,integrate_affy_associations,relationship_type,'over-write previous')
                    index+=1
    print 'Finished integrating updated WikiPathways'

def importKEGGAssociations(selected_species,force):
    supported_databases = ['Ag','At','Ce','Dm','Dr','Hs','Mm','Os','Rn']
    getSourceData()
    
    if selected_species != None: ### Restrict by selected species
        supported_databases2=[]
        for species in selected_species:
            if species in supported_databases:
                supported_databases2.append(species)
        supported_databases = supported_databases2

    for species in supported_databases:
        if force == 'yes':
            downloadKEGGPathways(species)
        for mod in mod_types:
            extractKEGGAssociations(species,mod,system_codes)

def importPathwayCommons(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support PathwayCommons update.' % string.join(original_species,',')
    else:
        if force == 'yes':
            downloadPathwayCommons()
    
        getSourceData()
    
        for species in selected_species:
            for mod in mod_types:
                extractGMTAssociations(species,mod,system_codes,'PathwayCommons')
    
def importTranscriptionTargetAssociations(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support Transcription Factor association update.' % string.join(original_species,',')
    else:
        ### No need to specify a species since the database will be added only to currently installed species
        if force == 'yes':
            source = downloadPAZARAssocations()
        if source == 'raw':
            x = importPAZARAssociations() ### builds the PAZAR TF-symbol associations from resource.csv files
        else:
            x = importPAZARcompiled() ### imports from pre-compiled/downloaded TF-symbol associations
            
        y = importAmandeusPredictions(force)
        z = dict(x.items() + y.items())
        exportSymbolRelationships(z,selected_species,'TFTargets','mapp')

def importPhenotypeOntologyData(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support Phenotype Ontology update.' % string.join(original_species,',')
    else:
        ### No need to specify a species since the database will be added only to currently installed species    
        if force == 'yes':
            downloadPhenotypeOntologyOBO()
            downloadPhenotypeOntologyGeneAssociations()
        x = importPhenotypeOntologyGeneAssociations()
        exportSymbolRelationships(x,selected_species,'MPhenoOntology','go')
    
def importDiseaseOntologyAssociations(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support Disease Ontology update.' % string.join(original_species,',')
    else:
        if force == 'yes':
            downloadDiseaseOntologyOBO()
            downloadDiseaseOntologyGeneAssociations(selected_species)
        x = importDiseaseOntologyGeneAssocations()
        exportSymbolRelationships(x,selected_species,'CTDOntology','go')

def importGOSlimAssociations(selected_species,force):
    if force == 'yes':
        downloadGOSlimOBO()
    transferGOSlimGeneAssociations(selected_species)
    
def importMiRAssociations(selected_species,force):
    supported_databases = unique.read_directory('/'+database_dir)
    if selected_species != None: ### Restrict by selected species
        supported_databases=selected_species

    missing_miR_associations=[]
    found_miR_associations=[]
    for species in supported_databases:
        if force == 'yes':
            try:
                fn = downloadMiRDatabases(species)
                found_miR_associations.append((species,fn))
            except Exception:
                missing_miR_associations.append(species)
                
    for (species,fn) in found_miR_associations:
        importMiRGeneAssociations(species,fn)
        
def importBioMarkerAssociations(selected_species,force):
    original_species = selected_species
    selected_species = considerOnlyMammalian(selected_species)
    if len(selected_species) == 0:
        print 'PLEASE NOTE: %s does not support BioMarker association update.' % string.join(original_species,',')
    else:
        if force == 'yes':
            downloadBioMarkers()
        x = importBioMarkerGeneAssociations()
        exportSymbolRelationships(x,selected_species,'BioMarkers','mapp')

def importDomainAssociations(selected_species,force):
    if force == 'yes':
        paths = downloadDomainAssociations(selected_species)
        for (species,path) in paths:
            path = string.replace(path,'.gz','.txt')
            importDomainGeneAssociations(species, path)

def considerOnlyMammalian(selected_species):
    supported_mammals = ['Am','Bt', 'Cf', 'Ch', 'Cj', 'Cp', 'Do', 'Ec', 'Ee', 'Et', 'Fc', 'Gg', 'Go', 'Hs',
                         'La', 'Ma', 'Md', 'Me', 'Mi', 'Ml', 'Mm', 'Oa', 'Oc','Og', 'Op', 'Pc', 'Pp',
                         'Pt', 'Pv', 'Rn', 'Sa', 'Ss', 'St', 'Tb', 'Tn', 'Tr', 'Ts', 'Tt', 'Vp']
    filtered_species=[]
    if selected_species == None:
        selected_species = unique.read_directory('/'+database_dir)
        
    for i in selected_species:
        if i in supported_mammals:
            filtered_species.append(i)
    return filtered_species

def buildInferrenceTables(selected_species):
    for species_code in selected_species:
        file_found = verifyFile(database_dir+'/'+species_code+'/uid-gene/Ensembl-Symbol'+'.txt') ### If file is present, the below is not needed
        if file_found == 'no':
            try: gene_associations.swapAndExportSystems(species_code,'Ensembl','EntrezGene') ### Allows for analysis of Ensembl IDs with EntrezGene based GO annotations (which can vary from Ensembl)
            except Exception: null=[] ### Occurs if EntrezGene not supported
    
            ### Build out these symbol association files
            try: gene_associations.importGeneData(species_code,('export','Ensembl'))
            except Exception: null=[] ### Occurs if EntrezGene not supported
            try: gene_associations.importGeneData(species_code,('export','EntrezGene'))
            except Exception: null=[] ### Occurs if EntrezGene not supported

def buildAccessoryPathwayDatabases(selected_species,additional_resources,force):
    global database_dir
    global program_type
    program_type,database_dir = unique.whatProgramIsThis()
    
    buildInferrenceTables(selected_species) ### Make sure these tables are present first!!!
        
    #print 'Attempting to update:', string.join(additional_resources,',')
    if 'KEGG' in additional_resources:
        try: importKEGGAssociations(selected_species,force)
        except Exception: print 'KEGG import failed (cause unknown)'
    if 'Transcription Factor Targets' in additional_resources:
        try: importTranscriptionTargetAssociations(selected_species,force)
        except Exception: print 'Transcription Factor Targets import failed (cause unknown)'
    if 'Phenotype Ontology' in additional_resources:
        try: importPhenotypeOntologyData(selected_species,force)
        except Exception: print 'Phenotype Ontology import failed (cause unknown)'
    if 'Disease Ontology' in additional_resources:
        try: importDiseaseOntologyAssociations(selected_species,force)
        except Exception: print 'Disease Ontology import failed (cause unknown)'
    if 'GOSlim' in additional_resources:
        try: importGOSlimAssociations(selected_species,force)
        except Exception: print 'GOSlim import failed (cause unknown)'
    if 'miRNA Targets' in additional_resources:
        try: importMiRAssociations(selected_species,force)
        except Exception: print 'miRNA Targets import failed (cause unknown)'
    if 'BioMarkers' in additional_resources:
        try: importBioMarkerAssociations(selected_species,force)
        except Exception: print 'BioMarkers import failed (cause unknown)'#,traceback.format_exc()
    if 'Domains' in additional_resources:
        try: importDomainAssociations(selected_species,force)
        except Exception: print 'Domains import failed (cause unknown)'
    if 'PathwayCommons' in additional_resources:
        try: importPathwayCommons(selected_species,force)
        except Exception: print 'PathwayCommons import failed (cause unknown)'
    if 'Latest WikiPathways' in additional_resources:
        try: importWikiPathways(selected_species,force)
        except Exception: print 'WikiPathways import failed (cause unknown)'
            
if __name__ == '__main__':
    selected_species = ['Mm']
    force = 'no'
    additional_resources=['Latest WikiPathways','PathwayCommons','Transcription Factor Targets','Domains','BioMarkers']
    additional_resources+=['miRNA Targets','GOSlim','Disease Ontology','Phenotype Ontology','KEGG']
    additional_resources=['Latest WikiPathways']
    buildAccessoryPathwayDatabases(selected_species,additional_resources,force)
    