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
import update; reload(update)
import OBO_import
import gene_associations

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
            
def downloadHMDBMetaboCardFlatFile():
    url = 'http://www.hmdb.ca/public/downloads/current/metabocards.zip'
    fln,status = update.download(url,'BuildDBs/HMDB/','')

def downloadKEGGPathwayIDs():
    url = 'ftp://ftp.genome.jp/pub/kegg/pathway/map_title.tab'
    fln,status = update.download(url,'BuildDBs/HMDB/','')
    
class HMBDInfo:
    def __init__(self,entry_data):
        self.kegg_list=[]; self.smpdb_list=[]; self.protein_names=[]
        for name in entry_data:
            if entry_data[name] == 'Not Available': entry_data[name] = ''
        self.hmdb_id = entry_data['hmdb_id']
        self.description = entry_data['description']
        self.name = entry_data['name']
        self.secondary_id = entry_data['secondary_id']
        self.iupac = entry_data['iupac']
        self.biocyc_id = entry_data['biocyc_id']
        self.cas_number = entry_data['cas_number']
        self.chebi_id = entry_data['chebi_id']
        self.pubchem_compound_id = entry_data['pubchem_compound_id']
        self.kegg_compound_id = entry_data['kegg_compound_id']
        
        for name in entry_data:
            if 'kegg_id' in name and 'pathway' in name:
                id = entry_data[name][3:]
                if id in kegg_pathways:
                    KEGG_pathway = kegg_pathways[id][0]
                    self.kegg_list.append(KEGG_pathway+':'+entry_data[name])
            elif 'smpdb' in name and 'pathway' in name:
                SMPDB_name = string.replace(name,'smpdb_id','name')
                SMPDB_pathway = entry_data[SMPDB_name]
                if len(entry_data[name])>0:
                    self.smpdb_list.append(SMPDB_pathway+':'+entry_data[name])
            elif 'metabolic_enzyme' in name and 'gene_name' in name:
                if len(entry_data[name])>0:
                    self.protein_names.append(entry_data[name])
        
    def HMDB(self): return self.hmdb_id
    def Description(self): return self.description
    def Name(self): return self.name
    def SecondaryIDs(self): return self.secondary_id
    def IUPAC(self): return self.iupac
    def CAS(self): return self.cas_number
    def CheBI(self): return self.chebi_id
    def KEGGCompoundID(self): return self.kegg_compound_id
    def KEGG(self): return self.kegg_list
    def SMPDB(self): return self.smpdb_list
    def Pathways(self): return self.kegg_list+self.smpdb_list
    def PathwaysStr(self): return string.join(self.Pathways(),',')
    def BioCyc(self):return self.biocyc_id
    def PubChem(self):return self.pubchem_compound_id
    def ProteinNames(self):return self.protein_names
    def ProteinNamesStr(self):return string.join(self.protein_names,',')
    def __repr__(self): return self.HMDB()
    
def importHMDBMetaboCardFlatFile():
    filename = 'BuildDBs/HMDB/metabocards.txt'
    fields_to_store = ['hmdb_id','description','name','secondary_id','iupac','biocyc_id','cas_number','chebi_id','kegg_compound_id','pubchem_compound_id','pathway_1_kegg_id']
    fields_to_store+= ['metabolic_enzyme_1_gene_name','metabolic_enzyme_1_swissprot_id','metabolic_enzyme_2_gene_name','metabolic_enzyme_2_swissprot_id']
    fields_to_store+= ['pathway_1_smpdb_id','pathway_2_smpdb_id','pathway_3_smpdb_id','pathway_1_name','pathway_2_name','pathway_3_name','pathway_2_kegg_id','pathway_3_kegg_id']
    fn=filepath(filename); field_data=''; field_name=''; entry_data={}; hmdb=[]; global kegg_pathways
    kegg_pathways = gene_associations.importGeneric('BuildDBs/HMDB/map_title.tab'); x=0
    for line in open(fn,'rU').xreadlines():
        data = cleanUpLine(line)
        if len(data)>0:
            if data[0]=='#': field_name = data[2:-1]
            else: field_data += data 
        else:
            if field_name in fields_to_store:
                entry_data[field_name] = field_data
            #else: print [field_name]
            field_name=''; field_data=''
        if 'END_METABOCARD' in data:
            ed = HMBDInfo(entry_data)
            hmdb.append(ed)
            entry_data={}
            x+=1
            #if x>5: break
    print len(hmdb),'HMDB entries obtained'
    exportTables(hmdb)

def exportTables(metabolite_list):
    infer_enzyme_to_metabolite_pathway_data = 'no'
    
    current_species_dirs = unique.returnDirectories('/Databases')
    ### Save results to all species directories
    for species_code in current_species_dirs:
        print 'Exporting metabolite data for:',species_code
        gene_dir = 'Databases/'+species_code+'/gene/HMDB.txt'
        gene_data = export.ExportFile(gene_dir)
        hmdb_cas_dir = 'Databases/'+species_code+'/uid-gene/HMDB-CAS.txt'
        hmdb_cas_data = export.ExportFile(hmdb_cas_dir)
        hmdb_chebi_dir = 'Databases/'+species_code+'/uid-gene/HMDB-ChEBI.txt'
        hmdb_chebi_data = export.ExportFile(hmdb_chebi_dir)
        hmdb_pubchem_dir = 'Databases/'+species_code+'/uid-gene/HMDB-PubChem.txt'
        hmdb_pubchem_data = export.ExportFile(hmdb_pubchem_dir)
        hmdb_keggcomp_dir = 'Databases/'+species_code+'/uid-gene/HMDB-KeggCompound.txt'
        hmdb_keggcomp_data = export.ExportFile(hmdb_keggcomp_dir)
        hmdb_mapp_dir = 'Databases/'+species_code+'/gene-mapp/HMDB-MAPP.txt'
        hmdb_mapp_data = export.ExportFile(hmdb_mapp_dir)
        cas_denom_dir = 'Databases/'+species_code+'/gene-mapp/denominator/CAS.txt'
        cas_denom_data = export.ExportFile(cas_denom_dir)
        hmdb_go_dir = 'Databases/'+species_code+'/gene-go/HMDB-GeneOntology.txt'
        if infer_enzyme_to_metabolite_pathway_data == 'yes':
            hmdb_go_data = export.ExportFile(hmdb_go_dir)
        
        headers = ['hmdb_id','name','description','secondary_id','iupac','cas_number','chebi_id','pubchem_compound_id','Pathways','ProteinNames']
        headers = string.join(headers,'\t')+'\n'
        gene_data.write(headers)
        
        ### Attempt to add GO and pathway data from database based on associated protein IDs (simple translation from human)
        mod = 'Ensembl'
        try: gene_annotations = gene_associations.importGeneData(species_code,mod)
        except Exception:
            mod = 'EntrezGene'
            try: gene_annotations = gene_associations.importGeneData(species_code,mod)
            except Exception: gene_annotations={}
        symbol_associations={}
        for geneid in gene_annotations: symbol_associations[gene_annotations[geneid].SymbolLower()] = geneid
        gotype = 'null'
        try: gene_to_go = gene_associations.importGeneGOData(species_code,mod,gotype)
        except Exception: gene_to_go={}
        try: gene_to_mapp = gene_associations.importGeneMAPPData(species_code,mod)
        except Exception: gene_to_mapp = {}
        
        for ed in metabolite_list:
            values = [ed.HMDB(),ed.Name(),ed.Description(),ed.SecondaryIDs(),ed.IUPAC(),ed.CAS(),ed.CheBI(),ed.PubChem(),ed.PathwaysStr(),ed.ProteinNamesStr()]
            values = string.join(values,'\t')+'\n'; gene_data.write(values)
            if len(ed.Pathways())>1:
                for pathway in ed.Pathways():
                    values = [ed.HMDB(),'Ch',pathway]; values = string.join(values,'\t')+'\n'; hmdb_mapp_data.write(values)
            
            if len(ed.CAS())>0:
                values = [ed.HMDB(),ed.CAS()]; values = string.join(values,'\t')+'\n'; hmdb_cas_data.write(values)
                values = [ed.CAS(),'Ca']; values = string.join(values,'\t')+'\n'; cas_denom_data.write(values)
            if len(ed.CheBI())>0: values = [ed.HMDB(),ed.CheBI()]; values = string.join(values,'\t')+'\n'; hmdb_chebi_data.write(values)
            if len(ed.PubChem())>0: values = [ed.HMDB(),ed.PubChem()]; values = string.join(values,'\t')+'\n'; hmdb_pubchem_data.write(values)
            if len(ed.KEGGCompoundID())>0: values = [ed.HMDB(),ed.KEGGCompoundID()]; values = string.join(values,'\t')+'\n'; hmdb_keggcomp_data.write(values)
            temp_go={}; temp_mapp={}
            
            if infer_enzyme_to_metabolite_pathway_data == 'yes':
                ### If associated enzyme annotated, use the gene symbol to find GO terms associated with the gene symbol for the metabolite
                ### Not sure if this is a bad idea or not
                for protein_name in ed.ProteinNames():
                    protein_name = string.lower(protein_name)
                    if protein_name in symbol_associations:
                        geneid = symbol_associations[protein_name]
                        if geneid in gene_to_go:
                            for goid in gene_to_go[geneid]: temp_go[goid]=[]
                        if geneid in gene_to_mapp:
                            for mapp in gene_to_mapp[geneid]: temp_mapp[mapp]=[]
                for goid in temp_go:
                    values = [ed.HMDB(),'GO:'+goid]; values = string.join(values,'\t')+'\n'; hmdb_go_data.write(values)     
                for mapp in temp_mapp:
                    values = [ed.HMDB(),'Ch',mapp]; values = string.join(values,'\t')+'\n'; hmdb_mapp_data.write(values)
                        
        gene_data.close(); hmdb_mapp_data.close(); hmdb_cas_data.close(); hmdb_chebi_data.close(); hmdb_pubchem_data.close();
        if infer_enzyme_to_metabolite_pathway_data == 'yes':
            hmdb_go_data.close()
        print 'File:',gene_dir,'exported.'
        
def verifyFile(filename):
    fn=filepath(filename); file_found = 'yes'
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception: file_found = 'no'
    return file_found
            
def buildMetabolomicsDatabase(force):
    ### No need to specify a species since the database will be added only to currently installed species
    if force == 'yes':
        downloadHMDBMetaboCardFlatFile()
        downloadKEGGPathwayIDs()
    importHMDBMetaboCardFlatFile()
    
if __name__ == '__main__':
    buildMetabolomicsDatabase('no'); sys.exit()
    downloadHMDBMetaboCardFlatFile()#;sys.exit()
    downloadKEGGPathwayIDs()
    importHMDBMetaboCardFlatFile()#;sys.exit()
    