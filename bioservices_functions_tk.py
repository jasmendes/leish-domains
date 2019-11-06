# -*- coding: cp1252 -*-

# https://datasciencelab.wordpress.com/2013/12/12/clustering-with-k-means-in-python/

# http://www.ncbi.nlm.nih.gov/protein/P43403

# http://www.ebi.ac.uk/Tools/webservices/tutorials/06_programming/python/soap/soappy

# http://pythonhosted.org/bioservices/references.html

#################################################################################	
##### Bioservices                                                        ########
#################################################################################
# https://pypi.python.org/pypi/bioservices
# http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan_soap
# Bioservices : a Python package to access biological web services programmatically¶

"""
# -*- python -*-
#
#  This file is part of bioservices software
#
#  Copyright (c) 2013-2014 - EBI-EMBL
#
#  File author(s):
#      Thomas Cokelaer <cokelaer@ebi.ac.uk>
#
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  website: https://github.com/cokelaer/bioservices
#  documentation: http://packages.python.org/bioservices
#
##############################################################################


Access to Biological Web Services from Python
Package Documentation
Bioservices is a Python package that provides access to many Bioinformatices Web Services (e.g., UniProt) and a framework to easily implement Web Services wrappers (based on WSDL/SOAP or REST protocols).
The primary goal of BioServices is to use Python as a glue language to provide a programmatic access to several Bioinformatics Web Services. By doing so, elaboration of new applications that combine several of the wrapped Web Services is fostered.
One of the main philosophy of BioServices is to make use of the existing biological databases (not to re-invent new databases) and to alleviates the needs for expertise in Web Services for the developers/users.
BioServices provides access to
###25 Web Services including:


BioModels
KEGG-
UniProt-
quickGO-
PSIQUIC
WikiPAthway
UniChem,
ChEMBL,
EUtils,
GeneProf,
PathwayCommons
BioDBNet,
UniChem).

Keywords:
BioServices,WebServices,Biology,BioDBNet,ChEBI,UniChem,Kegg,KEGG,BioModels,
EUtils,UniProt,PICR,ArrayExpress,MUSCLE,QuickGO,PDB,PSICQUIC,Blast,BioMART,
BioGRID,MIRIAM,BioMart,GeneProf,EUtils,ChEMBL,ChemSpider,HGNC,PathwayCommons
"""


#import urllib2

# file to be written to
#file = "downloaded_file.html"

#url = "http://www.pythonforbeginners.com/"
#response = urllib2.urlopen(url)

#open the file for writing
#fh = open(file, "w")

# read from request while writing to file
#fh.write(response.read())
#fh.close()

################################################################################
##### Loading libraries                                                    #####
################################################################################
from __future__ import division
import datetime, time
time_begin = datetime.datetime.fromtimestamp(time.time())

import sys, os
from glob import *
#print "FILE PATH:", sys.argv[0]
#print "FILE NAME:", os.path.basename(sys.argv[0])
#print sys.path
sys.path.append("C:\Python273\Lib\site-packages")
sys.path.append("C:\Python279\Lib\site-packages")
sys.path.append("C:\Python276\Lib\site-packages")
sys.path.append("C:\Python276\Lib")
#print sys.path




from bioservices import *

from Bio import Entrez
from Bio import ExPASy
from Bio import Medline
from Bio import PDB
from Bio import SwissProt
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbipsiblastCommandline
from Bio.Align.Applications import ClustalwCommandline

from Bio import AlignIO
from bioservices import QuickGO
from math import pi
import urllib


Entrez.email = 'fc34880@alunos.fc.ul.pt'
#Entrez.email = "jmend3z@gmail.com"

#import iprscan_soappy2
import iprscan_urllib2
#from  task07 import *

import random
import winsound
# ValueError: frequency must be in 37 thru 32767

import collections
from collections import Counter

import numbers
import decimal
import matplotlib as plt

from scipy import stats


global theInputFiles
theInputFiles = []



################################################################################

###############################################################################3
######   BioPython   -   TAX ID
###############################################################################

# http://stackoverflow.com/questions/16504238/attempting-to-obtain-taxonomic-information-from-biopython

import sys

def get_tax_id(species):
    
    """to get data from ncbi taxomomy, we need to have the taxid. we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(' ', "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    return record['IdList'][0]

def get_tax_data(taxid):
    
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)


def from_name_get_tax_id(species_list):

    print ("\nSearching Taxonomic name . . . ")
    #mail = raw_input("Enter email address : ")
    #while len(mail.strip().split("@")) != 2:
    #    mail = raw_input("Please Enter email address : ")

    #mail = Entrez.email
    
    #if not Entrez.email:
    #print "you must add your email address"
    #sys.exit(2)
    
    taxid_list = []  # Initiate the lists to store the data to be parsed in
    data_list = []
    lineage_list = []
    print('parsing taxonomic data...\n') # message declaring the parser has begun

    if species_list == []:
        print ("INFO:\nUse text box to Enter a list of tax names and click ready, then use tax name button. ")
        return 

    for species in species_list:
        print ('\t> '+species) # progress messages
        taxid = get_tax_id(species) # Apply your functions
        print ("Specie Name : ",species,"\tTax ID : ", taxid)
        data = get_tax_data(taxid)
        #print species, data
        lineage = {d['Rank']:d['ScientificName'] for d in data[0]['LineageEx'] if d['Rank'] in ['phylum']}
        print (species, lineage)

        reply = raw_input("Show more info ? (Y/N) : ")
        if reply == "Y":
            print (species, "\n",data)
        elif reply == "y":
            print (species, "\n",data)
        elif reply == "yes":
            print (species, "\n",data)
        else:
            print ("Thanks")
            
        taxid_list.append(taxid) # Append the data to lists already initiated
        data_list.append(data)
        lineage_list.append(lineage)

    print('complete!')


#species_list = ["Homo sapiens",'Helicobacter pylori 26695', 'Thermotoga maritima MSB8', 'Deinococcus radiodurans R1', 'Treponema pallidum subsp. pallidum str. Nichols', 'Aquifex aeolicus VF5', 'Archaeoglobus fulgidus DSM 4304']
#species_list = ['Helicobacter pylori', 'Thermotoga maritima ', 'Deinococcus radiodurans R1', 'Treponema pallidum subsp. pallidum str. Nichols', 'Aquifex aeolicus VF5', 'Archaeoglobus fulgidus DSM ']
#species_list = ['26695', 'Thermotoga maritima MSB8', 'Deinococcus radiodurans R1', 'Treponema pallidum subsp. pallidum str. Nichols', 'Aquifex aeolicus VF5', 'Archaeoglobus fulgidus DSM 4304']

#from_name_get_tax_id(species_list)


    
##########################################################################################
### BIOSERVICES FUNCTIONS
##########################################################################################

############################################################################################################
################################################################################
#####  UniProt                                                            ######
################################################################################
"""
The following status codes may be returned:

Code	Description
200	The request was processed successfully.
400	Bad request. There is a problem with your input.
404	Not found. The resource you requested doesn't exist.
410	Gone. The resource you requested was removed.
500	Internal server error. Most likely a temporary problem, but if the problem persists please contact us.
503	Service not available. The server is being updated, try again later."""
################################################################################



def from_tax_id_get_bs_uniprot(tax_id):

        time_begin = datetime.datetime.fromtimestamp(time.time())
        

        ### SEARCH BY TAX_ID 
        u = UniProt(verbose=False)
        # print dir(u)
        print ("""Database: Bioservices - UniProt
        Columns:
        Correct values are ['citation', 'clusters', 'comments', 'database', 'domains', 'domain', 'ec', 'id', 'entry name',
        'existence', 'families', 'features', 'genes', 'go', 'go-id', 'interpro', 'interactor', 'keywords', 'keyword-id',
        'last-modified', 'length', 'organism', 'organism-id', 'pathway', 'protein names', 'reviewed', 'score', 'sequence', '3d',
        'subcellular locations', 'taxonomy', 'tools', 'version', 'virus hosts']""")
        data = u.search("taxonomy:"+str(tax_id), frmt="tab", limit=1000000000, columns="id, entry name, protein names, interpro, go-id, go, genes,ec , pathway, subcellular locations, existence, organism, organism-id")
        res = data.split("\n")
        print ("\n[1] Taxonomy Filter \n", data.split("\t")[0], tax_id)
        
        for i in res[0:10]:
            print (i)
        dig  = (random.randint(0,1000))
            
        print ("Proteins found : ", len(res)-3)
        aFile = open("output\\tax"+str(tax_id)+"_"+str(dig)+"_bs_uniprot.txt", "wt") #"_ac_entry_ipr_name_gene_go_path_local_exi.txt", "wt")
        aFile.write(data)
        aFile.close()
        print (aFile, " saved!")
        theInputFiles.append(aFile.name)
        print (theInputFiles)
        print ("\n\tBioservices - UniProt Access TAX - DONE!\n")

        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print ("(complete!)")
        
        time_end = datetime.datetime.fromtimestamp(time.time())
        print("Time elapsed: ", str(time_end - time_begin))
        
        #return aFile
        
            
#from_tax_id_get_bs_uniprot("5671")


def from_uni_id_get_bs_uniprot1(uni_ids):
    
        ### SEARCH BY UNI_ID
        time_begin = datetime.datetime.fromtimestamp(time.time())

        from bioservices import UniProt
        u = UniProt(verbose=False)
        print ("""
        Database: Bioservices - UniProt
        Correct values are ['citation', 'clusters', 'comments', 'database', 'domains', 'domain', 'ec', 'id', 'entry name',
        'existence', 'families', 'features', 'genes', 'go', 'go-id', 'interpro', 'interactor', 'keywords', 'keyword-id',
        'last-modified', 'length', 'organism', 'organism-id', 'pathway', 'protein names', 'reviewed', 'score', 'sequence', '3d',
        'subcellular locations', 'taxonomy', 'tools', 'version', 'virus hosts']""")
    
        #s = ",".join(uni_ids)
        count = 0
        name = random.choice(uni_ids)
        aFile = open("output\\"+str(name)+"_"+str(len(uni_ids))+"_bs_uniprot.txt", "w") # _entry_ipr_name_gene_go_path_local_exi
        #aFile.write("id\t  protein names\t interpro\t go-id\t go\t genes\t pathway\t subcellular\t locations\t existence\t organism\t organism-id\n")
        for i in uni_ids:
            data = u.search(i, frmt="tab", limit=1000000000, columns="id,  protein names, interpro, domain, domains, go-id, go, genes,ec , pathway, subcellular locations, existence, organism, organism-id")#" organism-id, entry name, id, protein names, interpro,  genes, go-id, go,  pathway, subcellular locations, existence, organism") # limit seacrh (2march2015) =153
            #print help(data)
            res1 = data.split("\n")
            print (count, i)
            if count == 0:
                for i in res1[:-1]:
                    #print i
                    aFile.write(i+"\n")
                    count += 1
            else:
                for i in res1[1:-1]:
                    #print i
                    aFile.write(i+"\n")
                    count += 1

        #for i in res1:
            #print i

        
            
        print ("Proteins found : ", count-1)
        aFile.close()
        print (aFile, " saved!")
        theInputFiles.append(aFile.name)

        # WITH UNI ID GET TAX ID PROTEINS

        af = open(str(aFile.name), "r")
        

        for i in af.readlines()[1:]:
            i = i.strip()
            #print i.split("\t")[-1]
            items = i.split("\t")
            if len(items) > 5:
                new_tax_id = items[-1]
                #print new_tax_id
                break

        #thefilename = from_tax_id_get_bs_uniprot(new_tax_id)

        #match_uni_ids_from_list(uni_ids, thefilename.name)
        
        print ("\n\tBioservices - UniProt Access ACS - DONE!\n")
        
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print ("(complete!)")

        time_end = datetime.datetime.fromtimestamp(time.time())
        print("Time elapsed: ", str(time_end - time_begin))

        
        

#uni_ids = ["A4I5T4","A4I0N4","A4I9M5"]#"A4HVK7","E9AHP1","A4HZH7","A4I8Y5","E9AI05","A4HRR9","A4I5T4","A4I0N4","A4I9M5","A4HVK7","E9AHP1","A4HZH7","A4I8Y5","E9AI05","A4HRR9"]
#from_uni_id_get_bs_uniprot1(uni_ids)

################################################################################

def from_gene_id_get_bs_uniprot1(gene_ids):
    
        ### SEARCH BY UNI_ID
        time_begin = datetime.datetime.fromtimestamp(time.time())

        from bioservices import UniProt
        u = UniProt(verbose=False)
        print ^("""
        Database: Bioservices - UniProt
        Correct values are ['citation', 'clusters', 'comments', 'database', 'domains', 'domain', 'ec', 'id', 'entry name',
        'existence', 'families', 'features', 'genes', 'go', 'go-id', 'interpro', 'interactor', 'keywords', 'keyword-id',
        'last-modified', 'length', 'organism', 'organism-id', 'pathway', 'protein names', 'reviewed', 'score', 'sequence', '3d',
        'subcellular locations', 'taxonomy', 'tools', 'version', 'virus hosts']""")
    
        #s = ",".join(uni_ids)
        count = 0
        name = random.choice(gene_ids)
        aFile = open("output\\"+str(name)+"_"+str(len(gene_ids))+"_bs_uniprot.txt", "w") # _entry_ipr_name_gene_go_path_local_exi
        #aFile.write("id\t  protein names\t interpro\t go-id\t go\t genes\t pathway\t subcellular\t locations\t existence\t organism\t organism-id\n")
        for i in gene_ids:
            data = u.search(i, frmt="tab", limit=1000000000, columns="id,  protein names, interpro, domain, domains, go-id, go, genes,ec , pathway, subcellular locations, existence, organism, organism-id")#" organism-id, entry name, id, protein names, interpro,  genes, go-id, go,  pathway, subcellular locations, existence, organism") # limit seacrh (2march2015) =153
            #print help(data)
            res1 = data.split("\n")
            print (count, i)
            if count == 0:
                for i in res1[:-1]:
                    #print i
                    aFile.write(i+"\n")
                    count += 1
            else:
                for i in res1[1:-1]:
                    #print i
                    aFile.write(i+"\n")
                    count += 1

        #for i in res1:
            #print i

        
            
        print ("Proteins found : ", count-1)
        aFile.close()
        print (aFile, " saved!")
        theInputFiles.append(aFile.name)

        # WITH UNI ID GET TAX ID PROTEINS

        af = open(str(aFile.name), "r")
        

        for i in af.readlines()[1:]:
            i = i.strip()
            #print i.split("\t")[-1]
            items = i.split("\t")
            if len(items) > 5:
                new_tax_id = items[-1]
                #print new_tax_id
                break

        #thefilename = from_tax_id_get_bs_uniprot(new_tax_id)

        #match_uni_ids_from_list(uni_ids, thefilename.name)
        
        print ("\n\tBioservices - UniProt Access ACS - DONE!\n")
        
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print ("(complete!)")

        time_end = datetime.datetime.fromtimestamp(time.time())
        print("Time elapsed: ", str(time_end - time_begin))

#############################################################################################################

def search_protein_name_in_org_bs_uniprot():
        
    # Make a open place to search in uniprot words

    ### SEARCH TERM COMBINATION
    u = UniProt(verbose=False)
    """
        Database: Bioservices - UniProt
        Correct values are
        ['citation', 'clusters', 'comments', 'database', 'domains', 'domain', 'ec', 'id', 'entry name', 'existence',
        'families', 'features', 'genes', 'go', 'go-id', 'interpro', 'interactor', 'keywords', 'keyword-id', 'last-modified', 'length', 'organism',
        'organism-id', 'pathway', 'protein names', 'reviewed', 'score', 'sequence', '3d', 'subcellular locations', 'taxonomy', 'tools', 'version',
        'virus hosts']"""
    
    # print "\n\tUNIPROT  - Search Term\n"
    # print "Accession via entry name (e.g., ZAP70_HUMAN) is faster than by Entry (e.g., P43403)"
    # theTerms0 = "ec:1.2.1.18 reviewed:no AND human\n"
    theTerms0 = 'zap70+AND+organism:9606'
    print ("Example: ", theTerms0)
    theTerms = raw_input("Enter terms to find in organisms : ")
    data1 = u.search(str(theTerms), frmt="tab", limit=10000000, columns="id,entry name,protein names, genes,go-id,go,length,organism,organism-id")
    res1 = data1.split("\n")
    theTerms = theTerms.split("+")

    aFile = open("output\\"+str(theTerms[0])+"_ac_names_len_genes_go_organism.txt", "wt")
    count = 0
    for i in res1:
            print (i)
            aFile.write(i+"\n")
            count += 1
    print ("Proteins found : ", count-2)
    aFile.close()
    print (aFile, " saved!")
    theInputFiles.append(aFile.name)
    
    print ("\n\tBioservices - UniProt Access NAME - DONE!\n")
    ############################################################################

#search_protein_name_in_org_bs_uniprot()


def from_seq_get_bs_uniprot():
    # make a open place to search in uniprot words
    ###########################################################################
    ###########################################################################
    ### SEARCH TERM COMBINATION
    from bioservices import UniProt
    u = UniProt(verbose=False)
    #u.get_fasta_sequence("P43403")
    #sequence = u.searchUniProtId("P43403", "fasta")
    # u.timeout = 10   # some queries are long and requires much more time; default is 1000 seconds
    """
        Database: Bioservices - UniProt
        Correct values are ['citation', 'clusters', 'comments', 'database', 'domains', 'domain', 'ec', 'id', 'entry name', 'existence',
        'families', 'features', 'genes', 'go', 'go-id', 'interpro', 'interactor', 'keywords', 'keyword-id', 'last-modified', 'length', 'organism',
        'organism-id', 'pathway', 'protein names', 'reviewed', 'score', 'sequence', '3d', 'subcellular locations', 'taxonomy', 'tools', 'version',
        'virus hosts']"""
    theSeq0 = "MITIDGNGAVASVAFRTSEVIAIYPITPSSTMAEQADAWAGNGLKNVWGDTPRVVEMQSEAGAIATVHGALQTGALSTSFTSSQGLLLMIPTLYKLAGELTPFVLHVAARTVATHALSIFGDHSDVMAVRQTGCAMLCAANVQEAQDFALISQIATLKSRVPFIHFFDGFRTSHEINKIVPLADDTILDLMPQVEIDAHRARALNPEHPVIRGTSANPDTYFQSREATNPWYNAVYDHVEQAMNDFSAATGRQYQPFEYYGHPQAERVIILMGSAIGTCEEVVDELLTRGEKVGVLKVRLYRPFSAKHLLQALPGSVRSVAVLDRTKEPGAQAEPLYLDVMTALAEAFNNGERETLPRVIGGRYGLSSKEFGPDCVLAVFAELNAAKPKARFTVGIYDDVTNLSLPLPENTLPNSAKLEALFYGLGSDGSVSATKNNIKIIGNSTPWYAQGYFVYDSKKAGGLTVSHLRVSEQPIRSAYLISQADFVGCHQLQFIDKYQMAERLKPGGIFLLNTPYSADEVWSRLPQEVQAVLNQKKARFYVINAAKIARECGLAARINTVMQMAFFHLTQILPGDSALAELQGAIAKSYSSKGQDLVERNWQALALARESVEEVPLQPVNPHSANRPPVVSDAAPDFVKTVTAAMLAGLGDALPVSALPPDGTWPMGTTRWEKRNIAEEIPIWKEELCTQCNHCVAACPHSAIRAKVVPPEAMENAPASLHSLDVKSRDMRGQKYVLQVAPEDCTGCNLCVEVCPAKDRQNPEIKAINMMSRLEHVEEEKINYDFFLNLPEIDRSKLERIDIRTSQLITPLFEYSGACSGCGETPYIKLLTQLYGDRMLIANATGCSSIYGGNLPSTPYTTDANGRGPAWANSLFEDNAEFGLGFRLTVDQHRVRVLRLLDQFADKIPAELLTALKSDATPEVRREQVAALRQQLNDVAEAHELLRDADALVEKSIWLIGGDGWAYDIGFGGLDHVLSLTENVNILVLDTQCYSNTGGQASKATPLGAVTKFGEHGKRKARKDLGVSMMMYGHVYVAQISLGAQLNQTVKAIQEAEAYPGPSLIIAYSPCEEHGYDLALSHDQMRQLTATGFWPLYRFDPRRADEGKLPLALDSRPPSEAPEETLLHEQRFRRLNSQQPEVAEQLWKDAAADLQKRYDFLAQMAGKAEKSNTD"
    print ("Example : \n", theSeq0)
    theSeq = raw_input("\nEnter sequence : ")
    data1 = u.search(str(theSeq), frmt="tab", limit=100000000, columns="id,protein names, database,sequence")
    res1 = data1.split("\n")

    aFile = open("output\\seq_ac_names_len_genes_go_organism.txt", "wt")
    count = 0
    for i in res1:
            print (i)
            aFile.write(i+"\n")
            count += 1
    print ("Proteins found : ", count)
    aFile.close()
    print (aFile, " saved!")
    print ("\n\tBioservices - UniProt Access SEQ - DONE!\n")
    
    ############################################################################

#from_seq_get_bs_uniprot()


################################################################################
##### BIOMART QUERY InterPro                                          TAXID    #######
##################################################################################

filename1= "" #bmipr_taxid - action2
filename2= "" #bmipr_ipr2go - auto
#filename3= "" #taxid2ipr2go - action merged ->> raw_input("Enter Unique File Name")






###########################################################3
        
def file2list(theFilename):
    """
['A4HT31',
'IPR008974',
'TRAF-like',
'UniProt/TrEMBL',
'Leishmania infantum',
'5671',
'IPR008974',
'GO:0005515',
'protein binding',
'function']"""

    aProteins = {}
    aFile = open(theFilename, "rt")
    for aLine in aFile.readlines()[1:]:
        aRow = aLine.strip().split("\t")# FORMAT1 = ac_ipr_name_db_taxn_taxid_ipr_go_goid_root
        #aRow = aLine.strip().split("\t")# FORMAT2 = 5671_ac_entry_ipr_name_gene_go_goid_path_local_exi_tax
        #print aRow
        if len(aRow)== 10: # list index out of range, key error, 
            aProteins[str(aRow[0])].append({'ipr_name':aRow[2],'uniprot':aRow[0], 'interpro':aRow[1], 'tax_id':aRow[5], 'go_id':aRow[7],'go_term':aRow[8], 'go_root':aRow[9]})
        else:
            aProteins[str(aRow[0])] = [{'ipr_name':aRow[2],'uniprot':aRow[0],'interpro':aRow[1], 'tax_id':aRow[5], 'go_id':aRow[7],'go_term':aRow[8], 'go_root':aRow[9]}]
    aFile.close()
    #return aProteins
    
    #af = open(theFilename, "rt")
    acs = []
    iprs = []
    go_ids = []
    go_terms = []
    go_roots = []
    for i in aFile.readlines()[1:]:
        print (i)
        ac = i.strip.split("\t")[0]
        if ac not in acs:
            acs.append(ac)
            
        ipr = i.strip.split("\t")[1]
        if ipr not in iprs:
            iprs.append(ipr)
            
        go_id = i.strip.split("\t")[7]
        if go_id not in go_ids:
            go_ids.append(go_id)

        go_term = i.strip.split("\t")[8]
        if go_term not in iprs:
            go_terms.append(go_term)

        go_root = i.strip.split("\t")[9]
        go.roots.append(go_root)
        
        
    
###uni_ids = file2list("output/5671_bmq_ac_ipr_go.txt")
###from_uni_id_get_bs_ac2ipr2go(uni_ids)

# http://pythonhosted.org/bioservices/_modules/bioservices/biomart.html

"""REACTOME example::

        s.lookfor("reactome")
        s.datasets("REACTOME")
        ['interaction', 'complex', 'reaction', 'pathway']

        s.new_query()
        s.add_dataset_to_xml("pathway")
        s.add_filter_to_xml("species_selection", "Homo sapiens")
        s.add_attribute_to_xml("pathway_db_id")
        s.add_attribute_to_xml("_displayname")
        xmlq = s.biomartQuery.get_xml()
        res = s.query(xmlq)



to retrieve filters available for a dataset:

        :param str dataset: e.g. oanatinus_gene_ensembl

        ::

            >>> s.filters("uniprot").split("\n")[1].split("\t")
            >>> s.filters("pathway")["species_selection"]
            [Arabidopsis thaliana,Bos taurus,Caenorhabditis elegans,Canis familiaris,Danio
            rerio,Dictyostelium discoideum,Drosophila melanogaster,Escherichia coli,Gallus
            gallus,Homo sapiens,Mus musculus,Mycobacterium tuberculosis,Oryza
            sativa,Plasmodium falciparum,Rattus norvegicus,Saccharomyces
            cerevisiae,Schizosaccharomyces pombe,Staphylococcus aureus N315,Sus
            scrofa,Taeniopygia guttata ,Xenopus tropicalis]"""

###############################################
#### BIOMART                                ####
##########################################33


def from_tax_id_get_bs_interpro2go(tax_id):
    
        time_begin = datetime.datetime.fromtimestamp(time.time())
        s = BioMart()
        ret = s.registry() # to get information about existing services aka databases
        #dir(s)
        """['__class__', '__delattr__', '__dict__', '__doc__', '__format__',
        '__getattribute__', '__hash__', '__init__', '__module__', '__new__',
        '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__',
        '__str__', '__subclasshook__', '__weakref__', '_biomartQuery', '_databases',
        '_debugLevel', '_display_names', '_easyXMLConversion', '_fixing_encoding',
        '_fixing_unicode', '_get_databases', '_get_displayNames',
        '_get_easyXMLConversion', '_get_hosts', '_get_level', '_get_names', '_get_url',
        '_get_valid_attributes', '_hosts', '_init', '_names', '_set_easyXMLConversion',
        '_set_level', '_set_url', '_url', '_valid_attributes', '_xml_example',
        'add_attribute_to_xml', 'add_dataset_to_xml', 'add_filter_to_xml', 'attributes',
        'checkParam', 'configuration', 'create_attribute', 'create_filter', 'databases',
        'datasets', 'debugLevel', 'displayNames', 'easyXML', 'easyXMLConversion',
        'filters', 'getUserAgent', 'get_xml', 'hosts', 'last_response', 'logging',
        'lookfor', 'name', 'names', 'new_query', 'onWeb', 'pubmed', 'query', 'registry',
        'request', 'requestPost', 'response_codes', 'url', 'urlencode',
        'valid_attributes', 'version']"""
        s.names      # alias to list of valid service names from registry
        "unimart" in s.names
        s = BioMart(verbose=False)
        # s.lookfor("interpro")#Candidate
        # s.datasets("prod-intermart_1")
        print ("\nExpected time = 2 min \n") #  1 min 45 sec)
        #xmlq = s.biomartQuery.get_xml()
        #res = s.query(xmlq)
        """
        ATRIBUTOS - protein
        ['protein_length', 'method_id', 'pos_from', 'entry_ac',
        'supermatch_entry_name', 'crc64', 'entry_name', 'protein_tax_id',
        'match_status', 'supermatch_start_coordinate', 'method_name',
        'supermatch_stop_coordinate', 'protein_database_name', 'supermatch_entry_ac',
        'protein_name', 'fragment', 'taxonomy_full_name', 'taxonomy_scientific_name',
        'protein_accession', 'md5', 'entry_short_name', 'pos_to', 'entry_type',
        'match_score', 'supermatch_entry_short_name', 'supermatch_entry_type',
        'method_database_name']

        ATRIBUTOS - entry
        ['protein_ac' 'method_id' 'parent_entry_short_name' 'pos_from' 'abstract' 'found_in_entry_id'
        'protein_length' 'found_in_entry_type' 'crc64' 'entry_name' 'child_entry_name'
        'protein_tax_id' 'match_status' 'contained_entry_name' 'go_root_term'
        'go_id' 'method_name' 'contained_entry_type' 'protein_database_name'
        'found_in_entry_short_name' 'cross_reference_name' 'cross_reference_ac'
        'child_entry_short_name' 'protein_name' 'go_term_name' 'fragment'
        'contained_entry_id' 'contained_entry_short_name'
        'taxonomy_full_name' 'taxonomy_scientific_name' 'found_in_entry_name'
        'child_entry_type' 'md5''entry_short_name''pos_to''entry_type'
        'database_name'  'match_score' 'parent_entry_name' 'parent_entry_type'
        'child_entry_id' 'parent_entry_id' 'entry_id' 'method_database_name']

        ATRIBUTOS - uniparc
        print s.attributes("uniparc")

        {'match_status' 'entry_ac' 'method_id' 'upi' 
        'pos_to' 'pos_from'  'match_score' 'method_name' 'crc64' 'length' 
        'entry_short_name' 'entry_name' 'entry_type' 'method_database_name' }
        """
        """
        FILTERS
        ['match_status', 'protein_class_filter', 'protein_modelorg_filter',
        'protein_length_greater_than', 'protein_name', 'protein_tax_id_filter',
        'fragment', 'entry_type', 'protein_phylum_filter', 'entry_ac',
        'method_database_name', 'protein_length_less_than', 'crc64', 'method_name',
        'entry_name', 'protein_kingdom_filter', 'method_id', 'protein_database_name',
        'protein_accession', 'md5']
        """
        ###########################################################################################
        ### TASK 01 - NCBI ID: protein_tax_id_filter DATASET = "protein"               #######
        ###########################################################################################

        s.lookfor("interpro")                           #Candidate:
        s.datasets("prod-intermart_1")                  # MART name:
        # retrieve datasets available for this mart
        # ['protein', 'entry', 'uniparc']
        s.add_dataset_to_xml("protein")

        s.add_attribute_to_xml("protein_accession")
        #s.add_attribute_to_xml("protein_name")
        s.add_attribute_to_xml("supermatch_entry_ac") #total = 5283
        s.add_attribute_to_xml("supermatch_entry_name")
        #s.add_attribute_to_xml("entry_type")
        s.add_attribute_to_xml("protein_database_name")
        #s.add_attribute_to_xml("protein_tax_id")
        s.add_attribute_to_xml("taxonomy_scientific_name")
          
        #s.add_attribute_to_xml("go_id")
        

        ### s.add_filter_to_xml("protein_length_greater_than", 1000)
        s.add_filter_to_xml("protein_tax_id_filter", str(tax_id) )

        xml_query = s.get_xml()
        res = s.query(xml_query)

        res = res.split("\n") # TOTAL = 22262, 22288(6out)

        header ="protein accession\tsupermatch entry ac\tsupermatch entry name\tprotein database name\ttaxonomy scientific name\n"
        bottom =  "# Biomart Query Hits: %s \nTaxonomic ID : %s"  % (len(res), tax_id) # "%s=%s" % (k, v)
        print (header)
        for i in res[0:10]:
            print (i)
        aFile = open("output\\tax"+str(tax_id)+"_bmq_ac_ipr.txt", "wt")#YOUR FILENAME
        aFile.write(header)

        aLines = []
        for aLine in res:
           aFile.write(aLine+"\n")
           aLines.append(aLine)

        
        accs = []
        acs = []
        ipprs = []
        iprs = []
        ipr_names = []
        ippr_names = []
        for i in aLines:
            i= i.split("\t")
            ac = i[0]
            accs.append(ac)
            if ac not in acs:
                acs.append(ac)
            ipr = i[1] if len(i) >= 2  else "null"
            ipprs.append(ipr)
            if ipr not in iprs:
                iprs.append(ipr)

            ipr_name = i[2] if len(i) >= 3  else "null"
            ippr_names.append(ipr_name)
            if ipr_name not in ipr_names:
                ipr_names.append(ipr_name)

        from collections import Counter
        c = Counter(accs)
        print ("\nTOTAL ACS found : ", len(accs)) # 6597
        print ("UNIQUE ACS : ", len(c), "\n", c.most_common(10)) # 2088
        for i in dict(c.most_common(10)):
            print (i, c[i])

        d  = Counter(ipprs)
        print ("\nTOTAL IPRs found : ", len(ipprs)) # 6597
        print ("UNIQUE IPRs : ", len(d), "\n", d.most_common(10)) # 3606
        for i in dict(d.most_common(10)):
            print (i, d[i])

        e  = Counter(ippr_names)
        print ("\nTOTAL IPRs NAMES found : ", len(ippr_names)) # 6597
        print ("UNIQUE IPR NAMES : ", len(e), "\n", e.most_common(10)) # 3606
        for i in dict(e.most_common(10)):
            print (i, e[i])

        print (bottom)

        aFile.close()
        print ("BioServices - InterPro by TAX ID| Hits : " ,len(res))

        print (aFile, " saved!")
        #theInputFiles.append(aFile.name)
        
        #def access_interpro2go_with_biomart():
        ########################################################################
        ### TASK 02 - GO TERMS - BIOMART - INTERPRO - DATASET = ENTRY
        #######################################################################
        #s = BioMart()
        #ret = s.registry() # to get information about existing services aka databases
        #dir(s)
        """['__class__', '__delattr__', '__dict__', '__doc__', '__format__',
        '__getattribute__', '__hash__', '__init__', '__module__', '__new__',
        '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__',
        '__str__', '__subclasshook__', '__weakref__', '_biomartQuery', '_databases',
        '_debugLevel', '_display_names', '_easyXMLConversion', '_fixing_encoding',
        '_fixing_unicode', '_get_databases', '_get_displayNames',
        '_get_easyXMLConversion', '_get_hosts', '_get_level', '_get_names', '_get_url',
        '_get_valid_attributes', '_hosts', '_init', '_names', '_set_easyXMLConversion',
        '_set_level', '_set_url', '_url', '_valid_attributes', '_xml_example',
        'add_attribute_to_xml', 'add_dataset_to_xml', 'add_filter_to_xml', 'attributes',
        'checkParam', 'configuration', 'create_attribute', 'create_filter', 'databases',
        'datasets', 'debugLevel', 'displayNames', 'easyXML', 'easyXMLConversion',
        'filters', 'getUserAgent', 'get_xml', 'hosts', 'last_response', 'logging',
        'lookfor', 'name', 'names', 'new_query', 'onWeb', 'pubmed', 'query', 'registry',
        'request', 'requestPost', 'response_codes', 'url', 'urlencode',
        'valid_attributes', 'version']"""
        #s.names      # alias to list of valid service names from registry
        "unimart" in s.names
        #s = BioMart(verbose=False)
        # s.lookfor("interpro")#Candidate
        # s.datasets("prod-intermart_1")
        #print "shadow"
        #xmlq = s.biomartQuery.get_xml()
        #res = s.query(xmlq)
        """
        ATRIBUTOS - protein
        ['protein_length', 'method_id', 'pos_from', 'entry_ac',
        'supermatch_entry_name', 'crc64', 'entry_name', 'protein_tax_id',
        'match_status', 'supermatch_start_coordinate', 'method_name',
        'supermatch_stop_coordinate', 'protein_database_name', 'supermatch_entry_ac',
        'protein_name', 'fragment', 'taxonomy_full_name', 'taxonomy_scientific_name',
        'protein_accession', 'md5', 'entry_short_name', 'pos_to', 'entry_type',
        'match_score', 'supermatch_entry_short_name', 'supermatch_entry_type',
        'method_database_name']

        ATRIBUTOS - entry
        ['protein_ac' 'method_id' 'parent_entry_short_name' 'pos_from' 'abstract' 'found_in_entry_id'
        'protein_length' 'found_in_entry_type' 'crc64' 'entry_name' 'child_entry_name'
        'protein_tax_id' 'match_status' 'contained_entry_name' 'go_root_term'
        'go_id' 'method_name' 'contained_entry_type' 'protein_database_name'
        'found_in_entry_short_name' 'cross_reference_name' 'cross_reference_ac'
        'child_entry_short_name' 'protein_name' 'go_term_name' 'fragment'
        'contained_entry_id' 'contained_entry_short_name'
        'taxonomy_full_name' 'taxonomy_scientific_name' 'found_in_entry_name'
        'child_entry_type' 'md5''entry_short_name''pos_to''entry_type'
        'database_name'  'match_score' 'parent_entry_name' 'parent_entry_type'
        'child_entry_id' 'parent_entry_id' 'entry_id' 'method_database_name']

        ATRIBUTOS - uniparc
        print s.attributes("uniparc")

        {'match_status' 'entry_ac' 'method_id' 'upi' 
        'pos_to' 'pos_from'  'match_score' 'method_name' 'crc64' 'length' 
        'entry_short_name' 'entry_name' 'entry_type' 'method_database_name' }
        """
        """
        FILTERS
        ['match_status', 'protein_class_filter', 'protein_modelorg_filter',
        'protein_length_greater_than', 'protein_name', 'protein_tax_id_filter',
        'fragment', 'entry_type', 'protein_phylum_filter', 'entry_ac',
        'method_database_name', 'protein_length_less_than', 'crc64', 'method_name',
        'entry_name', 'protein_kingdom_filter', 'method_id', 'protein_database_name',
        'protein_accession', 'md5']
        """
        s.new_query()
        s.lookfor("interpro")                           #Candidate:
        s.datasets("prod-intermart_1")                  # MART name:)
        ### s.add_dataset_to_xml("protein")
        ### s.add_attribute_to_xml("protein_accession")
        s.add_dataset_to_xml("entry")
        #s.add_attribute_to_xml("protein_ac")
        s.add_attribute_to_xml("entry_id")#entry_name
        s.add_attribute_to_xml("entry_name")
        s.add_attribute_to_xml("entry_type")
        s.add_attribute_to_xml("go_id")
        s.add_attribute_to_xml("go_term_name")
        s.add_attribute_to_xml("go_root_term")
        #s.add_attribute_to_xml("database_name")
       
        #s.add_filter_to_xml("protein_tax_id_filter", str(tax_id) )
        xml_query = s.get_xml()
        res2 = s.query(xml_query)
        res2 = res2.split("\n") # TOTAL = 42128, 62706(6out)
        header22 ="""entry_id\tentry_name\tentry_type\tgo_id\tgo_term_name\tgo_root_name\n"""

        print (header22)
        for i in res2[0:10]: print (i)
        bFile = open("output\\bmq_ipr_go_"+str(len(res2))+".txt", "wt")
        bFile.write(str(header22))
        for aLine in res2:
                bFile.writelines(aLine+"\n")
        bFile.close()
        print  ("""Biomart Query | Hits : """, str(len(res2)))
        print (bFile, " saved!")
        
        #match_ac2_ipr2go(tax_id)
        #########
        #def match_ac2_ipr2go(tax_id):
        """('Biomart Query\n        Taxonomic ID """ +tax_id+ """  protein_accession    supermatch_entry_ac    supermatch_entry_name    protein_database_name    tax_id')
        O00913	IPR001951	Histone H4	UniProt/TrEMBL
        BioServices - InterPro | Hits :  XXXXX
        ('Biomart Query | Hits : 42128
        entry_id    entry_name    go_id    go_term_name    go_root_name  \n        ')
        IPR024898	Octanoyltransferase LipM	GO:0016415	octanoyltransferase activity	function"""

        af = open("output/tax"+str(tax_id)+"_bmq_ac_ipr.txt", "rt")
        bf = open("output/bmq_ipr_go_"+str(len(res2))+".txt", "rt")
        cf = open("output/tax"+str(tax_id)+"_bmq_ac_ipr_go.txt", "wt")# IPR are the keys
        #ef = open("output/"+str(tax_id)+"_bmq_ac_ipr.txt", "rt")
        #ef = open("output/"+str(tax_id)+"_bmq_ac_ipr.txt", "wt")
        

        
        header3 ="protein_accession\tsupermatch_entry_ac\tsupermatch_entry_name\tprotein_database_name\ttaxonomy_scientific_name\tentry_name\tentry_type\tgo_id\tgo_term_name\tgo_root_name\t\n"
        result = []
        count1 = 0
        count2 = 0
        acss = []
        iprss = []
        cf.write(header3)
        from collections import defaultdict
        ipr_dict1 = defaultdict(list)
        ipr_dict2 = defaultdict(list)
        for i1 in af.readlines()[2:]:
                i1 = i1.strip()
                #print i1
                result.append(i1+"||")
                ac1 = i1.split("\t")[0]
                if ac1 not in acss:
                    acss.append(ac1)
                #print ac1
                
                #del ipr[:]
                ipr1 = i1.split("\t")[1] if len(i1.split("\t")) > 1 else 'None'
                ipr_dict1[ipr1].append(i1) # print ac1
                #print ipr1
                if ipr1 not in iprss:
                    iprss.append(ipr1)
                
        for i2 in bf.readlines():
                i2 = i2.strip()
                #print i
                ipr2 = i2.split("\t")[0]
                ipr_dict2[ipr2] = i2[10:] # print ac1
                
        dd = defaultdict(list)
        for k,v in ipr_dict1.items():
            for i in v:
                #dd[k].append(i)
                if k in ipr_dict2.keys():
                        #print "igual" ,  "ipr2 = ipr1"
                        ##print i, ipr_dict1[i]
                        #print i, ipr_dict2[i]
                        count1 += 1
                        result.append(k)
                        merged = i +"\t"+ ipr_dict2[k]
                        merged = merged.replace("\n","\t")
                        #print merged
                        #print i1, "\n\t",i2
                        cf.write(merged+"\n")
                else:
                        #print ac1, ipr1, "diferent", ipr2 # 42k units to search each value
                        #print i1
                        merged = ipr_dict1[k]
                        #print merged
                        #merged = merged.replace("\n","\t")
                        cf.write(str(merged)+"\n")
                        
                        count2 += 1
                                
                #else:
                        #print "oi"
                        #print "Bioservices - InterPro - MERGE taxid2ipr2go - DONE!"

        print ("\nTOTAL Proteins : ", len(acss))
        print ("TOTAL IPR append :", count1)
        print ("TOTAL IPR not append : ", count2)
        print ("TOTAL IPRS in List : ", count1+count2)
        print ("Unique IPRS in List 1 and  2 : ", len(ipr_dict1.keys()),  len(ipr_dict2.keys()))
        print ("\n\tBioservices - InterPro - ACCESS AC2IPR2GO -  DONE!\n") # time = 1min 45 sec
        print (cf, " saved!")
        theInputFiles.append(cf.name)
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print ("(complete!)")
        
        time_end = datetime.datetime.fromtimestamp(time.time())
        print("Time elapsed: ", str(time_end - time_begin))

        #return cf.name

        
       
        

        
#match_ac2_ipr2go("5671")
#from_tax_id_get_bs_interpro2go("5671")
        



def from_uni_id_get_bs_ac2ipr2go(uni_ids):

        time_begin = datetime.datetime.fromtimestamp(time.time())
    
        s = BioMart()
        ret = s.registry() # to get information about existing services aka databases
        #dir(s)
        """['__class__', '__delattr__', '__dict__', '__doc__', '__format__',
        '__getattribute__', '__hash__', '__init__', '__module__', '__new__',
        '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__',
        '__str__', '__subclasshook__', '__weakref__', '_biomartQuery', '_databases',
        '_debugLevel', '_display_names', '_easyXMLConversion', '_fixing_encoding',
        '_fixing_unicode', '_get_databases', '_get_displayNames',
        '_get_easyXMLConversion', '_get_hosts', '_get_level', '_get_names', '_get_url',
        '_get_valid_attributes', '_hosts', '_init', '_names', '_set_easyXMLConversion',
        '_set_level', '_set_url', '_url', '_valid_attributes', '_xml_example',
        'add_attribute_to_xml', 'add_dataset_to_xml', 'add_filter_to_xml', 'attributes',
        'checkParam', 'configuration', 'create_attribute', 'create_filter', 'databases',
        'datasets', 'debugLevel', 'displayNames', 'easyXML', 'easyXMLConversion',
        'filters', 'getUserAgent', 'get_xml', 'hosts', 'last_response', 'logging',
        'lookfor', 'name', 'names', 'new_query', 'onWeb', 'pubmed', 'query', 'registry',
        'request', 'requestPost', 'response_codes', 'url', 'urlencode',
        'valid_attributes', 'version']"""
        s.names      # alias to list of valid service names from registry
        "unimart" in s.names
        s = BioMart(verbose=False)
        # s.lookfor("interpro")#Candidate
        # s.datasets("prod-intermart_1")
        print ("shadow")
        #xmlq = s.biomartQuery.get_xml()
        #res = s.query(xmlq)
        """
        ATRIBUTOS - protein
        ['protein_length', 'method_id', 'pos_from', 'entry_ac',
        'supermatch_entry_name', 'crc64', 'entry_name', 'protein_tax_id',
        'match_status', 'supermatch_start_coordinate', 'method_name',
        'supermatch_stop_coordinate', 'protein_database_name', 'supermatch_entry_ac',
        'protein_name', 'fragment', 'taxonomy_full_name', 'taxonomy_scientific_name',
        'protein_accession', 'md5', 'entry_short_name', 'pos_to', 'entry_type',
        'match_score', 'supermatch_entry_short_name', 'supermatch_entry_type',
        'method_database_name']
        

        ATRIBUTOS - entry
        ['protein_ac' 'method_id' 'parent_entry_short_name' 'pos_from' 'abstract' 'found_in_entry_id'
        'protein_length' 'found_in_entry_type' 'crc64' 'entry_name' 'child_entry_name'
        'protein_tax_id' 'match_status' 'contained_entry_name' 'go_root_term'
        'go_id' 'method_name' 'contained_entry_type' 'protein_database_name'
        'found_in_entry_short_name' 'cross_reference_name' 'cross_reference_ac'
        'child_entry_short_name' 'protein_name' 'go_term_name' 'fragment'
        'contained_entry_id' 'contained_entry_short_name'
        'taxonomy_full_name' 'taxonomy_scientific_name' 'found_in_entry_name'
        'child_entry_type' 'md5' 'entry_short_name''pos_to' 'entry_type'
        'database_name'  'match_score' 'parent_entry_name' 'parent_entry_type'
        'child_entry_id' 'parent_entry_id' 'entry_id' 'method_database_name']

        ATRIBUTOS - uniparc
        print s.attributes("uniparc")

        {'match_status' 'entry_ac' 'method_id' 'upi' 
        'pos_to' 'pos_from'  'match_score' 'method_name' 'crc64' 'length' 
        'entry_short_name' 'entry_name' 'entry_type' 'method_database_name' }
        """
        """
        FILTERS
        ['match_status', 'protein_class_filter', 'protein_modelorg_filter',
        'protein_length_greater_than', 'protein_name', 'protein_tax_id_filter',
        'fragment', 'entry_type', 'protein_phylum_filter', 'entry_ac',
        'method_database_name', 'protein_length_less_than', 'crc64', 'method_name',
        'entry_name', 'protein_kingdom_filter', 'method_id', 'protein_database_name',
        'protein_accession', 'md5']
        """
        ###########################################################################################
        ### TASK 01 - NCBI ID: protein_tax_id_filter DATASET = "entry"               #######
        ###########################################################################################

        s.lookfor("interpro")                           #Candidate:
        s.datasets("prod-intermart_1")                  # MART name:
        # retrieve datasets available for this mart
        # ['protein', 'entry', 'uniparc']
        s.add_dataset_to_xml("entry")

        #s._valid_attributes()
        #s.add_attribute_to_xml("protein_accession")
        s.add_attribute_to_xml("protein_ac")
        s.add_attribute_to_xml("protein_name")
        
        s.add_attribute_to_xml("entry_id") #total = 5283
        s.add_attribute_to_xml("entry_name")
        s.add_attribute_to_xml("entry_type")
        #s.add_attribute_to_xml("match_status") # 'entry_name
        s.add_attribute_to_xml("go_id")
        s.add_attribute_to_xml("go_term_name")
        s.add_attribute_to_xml("go_root_term")
        
        s.add_attribute_to_xml("protein_database_name")
        s.add_attribute_to_xml("taxonomy_scientific_name")
        s.add_attribute_to_xml("protein_tax_id")
        
        #s.add_attribute_to_xml("contained_entry_short_name")
        ### s.add_filter_to_xml("protein_length_greater_than", 1000)

        header = """protein_accession\tprotein_name\tentry_ac\tentry_name\tentry type\tgo_id\tgo_term\tgo_root\tprotein_database_name\ttax_name\ttax_id\n"""
        print (header)
        name = random.choice(uni_ids)
        aFile1 = open("output\\"+str(name)+"_"+str(len(uni_ids))+"_bs_biomart.txt", "wt") # +str(uni_ids[0])+"_bmq_ac_ipr_go_db.txt", "wt")#YOUR FILENAME
        aFile1.write(str(header))
        
        count1 = 0
        count2 = 0
        acs_no = []
        acs_yes = []
        for i in uni_ids: # Filter to keep it faster ~ 30 min
            #print i
            s.add_filter_to_xml('protein_ac', str(i) )
            xml_query = s.get_xml()
            res = s.query(xml_query)
            res = res.split("\n") # TOTAL = 22262, 22288(6out)
            
            #print res
            #aFile.write(res)
            #for i in res[0:10]: print i
            #res.pop()
            aLines = []
            for aLine in res:
                if len(res) == 1:
                    #print i, " has no match in InterPro-entry."
                    aFile1.write(i+"\n")
                    #count1 += 1
                    acs_no.append(i)
        
                else:
                    #res.pop()                    
                    #print aLine
                    if len(aLine) > 1:
                        if aLine not in aLines:
                            aFile1.write(aLine+"\n")
                            aLines.append(aLine)
                            #count2 += 1
                            if i not in acs_yes:
                                acs_yes.append(i)
                        
                    else:
                        continue
                    
        #res = res.split("\n") # TOTAL = 22262, 22288(6out)
        for i in res[0:10]: print (i)

        aFile1.close()
        print ("Total Ids for query : ", len(uni_ids)) # 16
        print ("BioServices - InterPro by UNI ID | match : " ,len(acs_yes)) # count2 # 34, 52, 43
        print ("BioServices - InterPro by UNI ID | no match : " , len(acs_no)) # count1 # 9
        print ("Coverage ", len(acs_no)/float(len(acs_yes)))

        print (aFile1, " saved!")

        # WITH UNI ID GET TAX ID PROTEINS

        af = open(str(aFile1.name), "r")
        

        for i in af.readlines()[1:]:
            i = i.strip()
            #print i.split("\t")[-1]
            items = i.split("\t")
            if len(items) > 5:
                new_tax_id = items[-1]
                #print new_tax_id
                break

        # thefilename = from_tax_id_get_bs_uniprot(new_tax_id)
        # match_uni_ids_from_list(uni_ids, thefilename.name)

        # MERGE 2 FILES
        
        #allf = from_tax_id_get_bs_interpro2go(new_tax_id)
        #selected_idsf = match_uni_ids_from_list(uni_ids, allf)
        #theInputFiles.append(selected_idsf)
        #theInputFiles.append(allf)
        

        print ("\n\tBioservices - InterPro - ACCESS AC2IPR2GO -  DONE!\n")
        
        theInputFiles.append(aFile1.name)
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print ("(complete!)")
        
        time_end = datetime.datetime.fromtimestamp(time.time())
        print("Time elapsed: ", str(time_end - time_begin))
        
        
#uni_ids = ["A4I9M5","A4HVK7","E9AHP1","A4HZH7","A4I8Y5","E9AI05","A4HRR9","A4I5T4","A4I0N4"] #,"A4I9M5","A4HVK7","E9AHP1","A4HZH7","A4I8Y5","E9AI05","A4HRR9"]

#from_uni_id_get_bs_ac2ipr2go(uni_ids)

        
def match_2_files():

    allf = from_tax_id_get_bs_interpro2go(new_tax_id)
    selected_idsf = match_uni_ids_from_list(uni_ids, allf)
    theInputFiles.append(selected_idsf.name)
    theInputFiles.append(allf.name)
        

##################################################################################################
######3 FROM UNI IDS GET TAX ID : THEN SEARCH LOCAL FILE FOR MATCH
##############################################################33

def from_uni_get_protein_tax(uni_ids):
    for i in uni_ids:
        print (i)

    



def from_ipr_id_get_bs_ipr2go(ipr_ids):
    
        ########################################################################
        ### TASK 02 - GO TERMS - BIOMART - INTERPRO - DATASET = ENTRY
        #######################################################################
        s = BioMart()
        ret = s.registry() # to get information about existing services aka databases
        #dir(s)
        """['__class__', '__delattr__', '__dict__', '__doc__', '__format__',
        '__getattribute__', '__hash__', '__init__', '__module__', '__new__',
        '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__',
        '__str__', '__subclasshook__', '__weakref__', '_biomartQuery', '_databases',
        '_debugLevel', '_display_names', '_easyXMLConversion', '_fixing_encoding',
        '_fixing_unicode', '_get_databases', '_get_displayNames',
        '_get_easyXMLConversion', '_get_hosts', '_get_level', '_get_names', '_get_url',
        '_get_valid_attributes', '_hosts', '_init', '_names', '_set_easyXMLConversion',
        '_set_level', '_set_url', '_url', '_valid_attributes', '_xml_example',
        'add_attribute_to_xml', 'add_dataset_to_xml', 'add_filter_to_xml', 'attributes',
        'checkParam', 'configuration', 'create_attribute', 'create_filter', 'databases',
        'datasets', 'debugLevel', 'displayNames', 'easyXML', 'easyXMLConversion',
        'filters', 'getUserAgent', 'get_xml', 'hosts', 'last_response', 'logging',
        'lookfor', 'name', 'names', 'new_query', 'onWeb', 'pubmed', 'query', 'registry',
        'request', 'requestPost', 'response_codes', 'url', 'urlencode',
        'valid_attributes', 'version']"""
        #s.names      # alias to list of valid service names from registry
        "unimart" in s.names
        #s = BioMart(verbose=False)
        # s.lookfor("interpro")#Candidate
        # s.datasets("prod-intermart_1")
        #print "shadow"
        #xmlq = s.biomartQuery.get_xml()
        #res = s.query(xmlq)
        """
        ATRIBUTOS - protein
        ['protein_length', 'method_id', 'pos_from', 'entry_ac',
        'supermatch_entry_name', 'crc64', 'entry_name', 'protein_tax_id',
        'match_status', 'supermatch_start_coordinate', 'method_name',
        'supermatch_stop_coordinate', 'protein_database_name', 'supermatch_entry_ac',
        'protein_name', 'fragment', 'taxonomy_full_name', 'taxonomy_scientific_name',
        'protein_accession', 'md5', 'entry_short_name', 'pos_to', 'entry_type',
        'match_score', 'supermatch_entry_short_name', 'supermatch_entry_type',
        'method_database_name']

        ATRIBUTOS - entry
        ['protein_ac' 'method_id' 'parent_entry_short_name' 'pos_from' 'abstract' 'found_in_entry_id'
        'protein_length' 'found_in_entry_type' 'crc64' 'entry_name' 'child_entry_name'
        'protein_tax_id' 'match_status' 'contained_entry_name' 'go_root_term'
        'go_id' 'method_name' 'contained_entry_type' 'protein_database_name'
        'found_in_entry_short_name' 'cross_reference_name' 'cross_reference_ac'
        'child_entry_short_name' 'protein_name' 'go_term_name' 'fragment'
        'contained_entry_id' 'contained_entry_short_name'
        'taxonomy_full_name' 'taxonomy_scientific_name' 'found_in_entry_name'
        'child_entry_type' 'md5''entry_short_name''pos_to''entry_type'
        'database_name'  'match_score' 'parent_entry_name' 'parent_entry_type'
        'child_entry_id' 'parent_entry_id' 'entry_id' 'method_database_name']

        ATRIBUTOS - uniparc
        print s.attributes("uniparc")

        {'match_status' 'entry_ac' 'method_id' 'upi' 
        'pos_to' 'pos_from'  'match_score' 'method_name' 'crc64' 'length' 
        'entry_short_name' 'entry_name' 'entry_type' 'method_database_name' }
        """
        """
        FILTERS
        ['match_status', 'protein_class_filter', 'protein_modelorg_filter',
        'protein_length_greater_than', 'protein_name', 'protein_tax_id_filter',
        'fragment', 'entry_type', 'protein_phylum_filter', 'entry_ac',
        'method_database_name', 'protein_length_less_than', 'crc64', 'method_name',
        'entry_name', 'protein_kingdom_filter', 'method_id', 'protein_database_name',
        'protein_accession', 'md5']
        """        
        # s.new_query()
        s.lookfor("interpro")                           #Candidate:
        s.datasets("prod-intermart_1")                  # MART name:)
        ### s.add_dataset_to_xml("protein")
        ### s.add_attribute_to_xml("protein_accession")
        s.add_dataset_to_xml("entry")
        #s.add_attribute_to_xml("protein_ac")
        s.add_attribute_to_xml("entry_id")#entry_name
        s.add_attribute_to_xml("entry_name")
        s.add_attribute_to_xml("entry_type")
        s.add_attribute_to_xml("go_id")
        s.add_attribute_to_xml("go_term_name")
        s.add_attribute_to_xml("go_root_term")
        #s.add_attribute_to_xml("protein_tax_id")
        
        #s.add_filter_to_xml("protein_tax_id_filter", str(tax_id) )
        xml_query = s.get_xml()
        res2 = s.query(xml_query)
        res2 = res2.split("\n") # TOTAL = 42128, 62706(6out)
        header2 = """entry_id\tentry_name\tentry type\tgo_id\tgo_term_name\tgo_root_name\n"""
        print (header2)
        for i in res2[0:50]: print (i)
        bFile = open("_bmq_ipr_go.txt", "wt")
        bFile.write(str(header2))
        for aLine in res2:
                bFile.writelines(aLine+"\n")
        bFile.close()
        print ("""Biomart Query | Hits : """ , str(len(res2)))
        print (bFile, " saved!\n")
        #########
        """('Biomart Query\n        InterPro ID List \n
        BioServices - InterPro | Hits :  22288
        ('Biomart Query | Hits : 42128\n       entry_id    entry_name    go_id    go_term_name    go_root_name  \n        ')
        IPR024898	Octanoyltransferase LipM	GO:0016415	octanoyltransferase activity	function"""
        #af = open("_bmq_ac_ipr.txt", "rt")
        bf = open("_bmq_ipr_go.txt", "rt")
        cf = open("_bmq_ipr_ipr2go.txt", "wt")# IPR are the keys
        header3 = "entry_id\tentry_name\tenrty type\tgo_id\tgo_term_name\tgo_root_name\n"
        cf.write(header3)
        result = []
        count1 = 0
        count2 = 0
        dict_iprs = {}
        from collections import defaultdict
        #dict_iprs = defaultdict(list)
        #ipr_ids= ipr_ids.split("\n")
        # ipr_ids= ipr_ids.split(",")
        # ipr_ids= ipr_ids.split(";")
        
        for i2 in bf.readlines():  
                i2 = i2.strip()
                #print i
                ipr2 = i2.split("\t")[0]
                rest = i2.split("\t")[1:]
                info = "\t".join(rest)
                dict_iprs[ipr2] = info
                #print rest
        #print dict_iprs.items()
        
                
        for i in ipr_ids:
                if i in dict_iprs.keys():
                        #print "igual" ,  "ipr2 = ipr1"
                        count1 += 1
                        # d[value].append(value)
                        #result=dict_iprs[i]
                        #"\t".join(result)
                        print (i,dict_iprs[i])
                        cf.write(i+"\t"+str(dict_iprs[i])+"\n")
                        #del result [:]
                else:
                        #print ac1, ipr1, "diferent", ipr2 # 42k units to search each value
                        #print i1
                        count2 += 1
                        #print "hello"
        print  (count1 , "/", count2, count1 + count2)
        print (cf, " saved!")
        theInputFiles.append(cf.name)

        print ("Bioservices  InterPro2GO - MATCH IPR2GO - DONE!")
        
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print ("(complete!)")
       

#ipr_ids = ["IPR024897", "IPR020605"]
#from_ipr_id_get_bs_ipr2go(ipr_ids)



        
#####################################################################################################################################
#####    QuickGO - Gene Ontology                                          ######
################################################################################

# http://pydoc.net/Python/bioservices/1.0.2/bioservices.quickgo/

def from_tax_id_get_quickgo(tax_id):
        
        aQG = QuickGO(verbose=False)
        # print dir(aQG)
        """['proteinDB', 'proteinID', 'proteinSymbol', 'qualifier', 'goID', 'goName', 'aspect', 'evidence', 'ref', 'with', 'proteinTaxon',
        'date', 'from', 'splice', 'proteinName', 'proteinSynonym', 'proteinType', 'proteinTaxonName', 'originalTermID', 'originalGOName']
        """
        aGoList = []
        #name = raw_input("\nSave filename as : ")
        aFile = open(str(tax_id)+"_uni_name_goid_type_db.txt","wt")
        count = 0
        aQueryResult = aQG.Annotation(tax=str(tax_id),col="proteinID,proteinName, proteinType, goID,goName,proteinSynonym, proteinDB,evidence, proteinTaxonName,proteinTaxon", frmt='tsv')
        for aLine in aQueryResult.splitlines():
            #print aLine
            acid1= aLine.split('\t')[0]
            #print str(id), aLine
            aGoList.append(str(aLine))
            aFile.write( str(aLine)+"\n")
            count += 1
        aFile.close()
        for i in aGoList[0:10]:
            print (i)

        print ("Proteins found for this TAX ",tax_id," : ", count)
        print ("\nBioservices - QuickGO - ACCESS TAX - DONE!")
        print (aFile , "is saved!")
        theInputFiles.append(aFile.name)
        
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print ("(complete!)")
        

#ids = ["5671"," "]
#for i in ids:
 #   from_tax_id_get_quickgo(i)

#####################################

def from_uni_ids_get_quickgo(uni_ids):
        
        aQG = QuickGO(verbose=False)
        # print dir(aQG)
        """['proteinDB', 'proteinID', 'proteinSymbol', 'qualifier', 'goID', 'goName', 'aspect', 'evidence', 'ref', 'with', 'proteinTaxon',
        'date', 'from', 'splice', 'proteinName', 'proteinSynonym', 'proteinType', 'proteinTaxonName', 'originalTermID', 'originalGOName']
        """
        aGoList = []
        name = uni_ids[0] # raw_input("\nSave filename as : ")
        aFile = open("output\\"+str(name)+"_"+str(len(uni_ids))+"_uni_goid_type_db.txt","wt")
        count = 0
        aLines = []
        for id in uni_ids:
                aQueryResult = aQG.Annotation(protein=id,col="proteinID,proteinName, proteinType, goID,goName,proteinTaxon,proteinTaxonName, proteinDB", frmt='tsv')
                # aQueryResult = aQG.Annotation(tax="5671",col="goID,goName", frmt='tsv')
                for aLine in aQueryResult.splitlines():
                        #print aLine
                        acid1= aLine.split('\t')[0]
                        #print str(id), aLine
                        aGoList.append(str(id)+ str(aLine))
                        count += 1
                        if aLine not in aLines:
                            #aLines.append(aLine)
                            aFile.write( str(aLine)+"\n")
                            
        aFile.close()
        print ("\nBioservices - QuickGO - ACCESS ACS - DONE!")
        print ("Proteins found for this AC List:", count)
        print (aFile , "is saved!")
        theInputFiles.append(aFile.name)

        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print ("(complete!)")
        

# ids = ["P12345","P12346"]
# from_uni_ids_get_quickgo(ids)
########################################################
                                                                                                # http://thomas-cokelaer.info/blog/2014/05/accessing-uniprot-with-python/

def from_go_ids_get_quickgo(go_ids):
        
        aQG = QuickGO(verbose=False)
        # print dir(aQG)
        """['proteinDB', 'proteinID', 'proteinSymbol', 'qualifier', 'goID', 'goName', 'aspect', 'evidence', 'ref', 'with', 'proteinTaxon',
        'date', 'from', 'splice', 'proteinName', 'proteinSynonym', 'proteinType', 'proteinTaxonName', 'originalTermID', 'originalGOName']
        """
        aGoList = []
        #name = raw_input("\nSave filename as : ")
        nam = go_ids[0]
        nam = nam.split(":")[-1]
        aFile = open(str(nam)+"_uni_goid_type_db.txt","wt")
        count = 0
        for id in go_ids:
                #aQueryResult = aQG.Annotation(protein=id,col="proteinID,proteinName, proteinType, goID,goName,proteinName,proteinTaxon,proteinTaxonName, proteinDB", frmt='tsv')
                aQueryResult = aQG.Annotation(goid=str(id),col="proteinID,proteinName, proteinType, goID,goName,proteinTaxon,proteinTaxonName, proteinDB", frmt='tsv')
                for aLine in aQueryResult.splitlines():
                        # print aLine
                        acid1= aLine.split('\t')[0]
                        # print str(id), aLine
                        aGoList.append(str(id)+ str(aLine))
                        aFile.write( str(aLine)+"\n")
                        count += 1
        aFile.close()
        print ("\nBioservices - QuickGO - ACCESS GOS - DONE!"        )
        print ("Proteins found for this GO term:", count)
        print (aFile , " is saved!")
        
        theInputFiles.append(aFile.name)
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print ("(complete!)")

# ids = ["GO:0008270"," "]
# from_go_ids_get_quickgo(ids)


def from_go_ids_get_quickgo_in_specie(go_ids):

    import re, sys
    from bioservices import QuickGO

    s = QuickGO(verbose=False)

    tax_id = raw_input("Enter Tax ID : ")
    # récupére les annotations pour l'apoptose
    res = s.Annotation(goid="GO:0006915", frmt="tsv", \
                       col="proteinDB,proteinID,goID,goName,proteinTaxon,proteinTaxonName", tax=str(tax_id))
    # découpage dans un array
    res = re.split("\n", res)

    print (len(res))
    # 10002
    # le premier élément contient l'entête
    # le dernier élément est vide

    proteins = []
    for r in res:
        cols = re.split("\t", r)
        # afficher les lignes pour UniProtKB
        if re.match('UniProtKB', cols[0]):
            proteins.append(cols[1])
            print (r)

    print (len(proteins))
    # 10000 protéines UniProtKB

    #sys.exit(0)
    print ("BIOSERVICES - Quick GO DONE!")
    theInputFiles.append(aFile.name)
    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    print ("(complete!)")

#ids = ["GO:0006915"]
#from_go_ids_get_quickgo_in_specie(ids)
    




#####################################################################################################################################
##### KEGG - KYOTO ENCYCLOPEDIA of GENES and GENOMES                      ####################################################3
#################################################################################################################################



def from_uni_ids_get_kegg_paths(uni_ids):

    
    from bioservices import KEGG
    k = KEGG(verbose=False)
    # ORGANISM : hsa
    # KEGG ID : hsa1234
    # PATHWAY : hsa12345    

    print ("\n\tKEGG  -  Look for relevant Pathways, using ac. . . ")

    from bioservices.uniprot import UniProt
    u = UniProt(verbose=False)
    kegg_ids = []
    organism_ids = []
    #filename = raw_input("Enter file name : ")
    namex = random.choice(uni_ids)
    co = len(uni_ids)
    af = open("kegg_data/"+str(namex)+"_"+str(co)+"_get_kegg_paths.txt","wt") # +str(filename)+
    af.write("UniProt ID \tKegg Organism \tKegg ID\tProtein Name \tKegg Pathways\n") #str(uni_i77d)+"\t"+str(k_organism)+"\t"+str(get_pname[:-2])+"\t"+str(kegg_paths1)+"\t"+str(paths))
    count = 1
    
    one_wait = 11 # seconds 2 = 27, 4=55, 6=1.8~107, 12=3.19~200, 179 =24:46
    expected_wait = len(uni_ids)*one_wait
    x = (expected_wait)/60
    g = float("{0:.1f}".format(x))
    print ("\nTotal ids: ",len(uni_ids)," >>>  Expected job time  : ", g, " minutes. . . ")
    for uni_id in uni_ids:
        
        print ("\n  Mapping ACC to KEGG ID . . . ")
        map1  = u.mapping(fr="ACC+ID", to="KEGG_ID", query=str(uni_id)) # AC
        #print len(map1), map1 # 1 defaultdict(<type 'list'>, {u'P43403': [u'hsa:7535']})
        for genes in map1.values():
            #i = str(i) # [u'lif:LINJ_10_0520'] or [u'lif:LINJ_28_2050', u'lif:LINJ_28_2060']
            #print genes # list
            for i in genes:
                #print  i
                kegg_id = i # # lif:LINJ_10_0520
                #print kegg_id
                kegg_ids.append(kegg_id)
                _id = i[4:] # LINJ_10_0520
                k_organism = i[0:3] # 'lif'
                #print uni_id, " > ", kegg_id #  > lif:LINJ_10_0520
                #kegg_ids.append(kegg_id)
        if k_organism not in organism_ids:
            organism_ids.append(k_organism)

        #print k.get(str(kegg_id))
        get_kegg = k.get(str(kegg_id)) # ENTRY (SEQ...)

        get_pname = k.find(str(k_organism), str(kegg_id))
        get_pname = str(get_pname) # NAMES

        print (count , uni_id, " > ",get_pname[:-1]) # lif:LINJ_14_0700	putative fatty acid elongase (EC:2.3.1.119)
        
            

        #print "Find pathway from gene . . . "
        #print help(map1)
        #print map1.values()
        #i = str(i)
        # print i

        get_path = k.get_pathway_by_gene(str(_id), str(k_organism)) # PATH
        # if no path is found, will print "No pathway found ?"
        
        #print "! No pathway found ! "
        #print get_path # if empty, will print None
        #paths = str(get_path)
        #dict2 = eval(str1)
        #print dict1==dict2
        #print type(get_path)
            
        
        if get_path.__repr__() == "None":
            #print get_path
            kegg_path = "None"
            af.write(str(uni_id)+"\t"+str(k_organism)+"\t"+str(get_pname[:-1])+"\t"+str(kegg_path)+"\n")
            
            

        elif isinstance(get_path, unicode):# isinstance('hello', str)
            #get_path ={}
            pathways = []
            
            kegg_path = get_path[0:8]
            print kegg_path
            kegg_path_name = get_path[10:]
            print kegg_path_name
            kegg_name_fix = kegg_path_name.replace("/","-")
            

            pathways.append(kegg_path+":"+kegg_path_name)
            # cf = open("kegg_data/"+str(kegg_path)+"_"+str(uni_id)+"_"+str(kegg_name_fix)+"_path.png","wb") # kegg_path_name
            # cf = open("kegg_data\\"+str(uni_id)+"_"+str(kegg_path)+"_"+str(kegg_name_fix)+"_path.png","wb") # 
            namex = "kegg_data\\"+str(kegg_path)+"_"+str(kegg_name_fix)+"_pathway.png"
            if os.path.isfile(namex) == False:
                cf = open(namex,"wb")
                res1 = []
                try:
                    # A4IAZ4
                    print "\nSaving image of the pathway  ",kegg_path, get_path[i]
                    #show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
                    res1 =  k.get(str(kegg_path)+"/image") # same as : res =  s.get("hsa05130","image")
                    cf.write(res1) # IMAGE
                except:
                    cf.close()

            kegg_paths = " |".join(pathways)
            print "\n> Pathway found: ", kegg_paths
            line = (str(uni_id)+"\t"+str(k_organism)+"\t"+str(get_pname[:-1])+"\t"+str(kegg_paths))
            line.replace("\n","\t")
            af.write(line+"\n")
                
            

        else:
            items = []
            pathways = []
            kegg_paths = []
            for i in get_path.values():
                #print i
                items.append(i)
            words = " |".join(items)
            #bf.write(str(words))
            for i in get_path.keys():
                if i.isdigit() == True and len(i) == 5:
                    #pathways.append(str(i)+":"+str(get_path[i]))# = str(i)
                    #print k_organism, i, map1[i]
                    kegg_path = i
                    kegg_path = k_organism + kegg_path
                    kegg_name_fix1 = get_path[i]
                    kegg_name_fix = kegg_name_fix1.replace("/","-")
                    
                    pathways.append(kegg_path+":"+get_path[i])
                    cf = open("kegg_data/"+str(uni_id)+"_"+str(kegg_path)+"_"+str(kegg_name_fix)+"_path.png","wb") # 
                    res1 = []
                    try:
                        # A4IAZ4
                        print "\nSaving image of the pathway  ",kegg_path, get_path[i]
                        #show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
                        res1 =  k.get(str(kegg_path)+"/image") # same as : res =  s.get("hsa05130","image")
                        cf.write(res1) # IMAGE
                    except:
                        cf.close()
            kegg_paths = " |".join(pathways)
            print "\n> Pathways found: ", kegg_paths
            #line = (str(uni_id)+"\t"+str(k_organism)+"\t"+str(get_pname[:-1])+"\t"+str(kegg_paths))
            #line.replace("\n","\t")
            
            #print kegg_paths
            ########################################   WRITE FILE   ###############################33
            allpaths = kegg_paths.split("|")
            for i in allpaths:
                line = (str(uni_id)+"\t"+str(k_organism)+"\t"+str(get_pname[:-1])+"\t"+str(i))
                af.write(line+"\n")
            

            # res0 =  k.get("hsa05130/image") # NF-kappa B signaling pathway - Homo sapiens (human)
            #print "Showing Pathway "# Marker in internet page . . ."
            
            #show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
                

        #print k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
        bf = open("kegg_data/"+str(uni_id)+"_kegg_entry.txt","w")
        
        bf.write(get_pname+"\n") # NAME
        bf.write(get_kegg) # ENTRY

        count +=1
        

    af.close()
    bf.close()
    print "Pathways found saved in" , af.name
    print "Pathways Info saved in:", bf.name 
    
    for i in organism_ids:
        df = open("kegg_data/"+str(i)+"_kegg_pathways.txt","wt")
        #print k.list("pathway", organism= i)
        res_paths = k.list("pathway", organism= i)
        pathways = [x.split()[0] for x in res_paths.strip().split("\n")]
        print len(pathways) ,"pathways found in Kegg organism ", i # as of Nov 2014 -> 94
        df.write(res_paths)
        df.close()
        print df, "saved!"
        
    

    #map2 = k.mapping(fr='KEGG_ID', to='ACC', format='tab', query='hsa:7535') # to="PDB_ID"
    #for i in map2.values():
    #print i#'P43403'
    print "Total ACS : ", count-1
    print "\nBioservices - Kegg - ACCESS ACS - DONE!"
    #print f1, "\n",f2,"\n", f3,"\n", f4 ,"\n", " are saved!"
    theInputFiles.append(af.name)
    
    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    print "(complete!)"

#onefile = open("2peps_250.txt","rt")

#ids = ["A4I0R1","A4HYI9"]
#from_uni_ids_get_kegg_paths(ids)

                            
    
#ids = ["A4IAZ4","A4HVQ0","A4I2L4"]#"P43403","E9AHJ2","A4HUV2","E9AG23","A4I3R9"]#"A4HXU4","A4HU65","Q9NGA0","A4HT41","A4HWD2","A4HWE9","E9AHJ2","A4I4C5"]
#ids = []

#for aLine in onefile.readlines():
 #       aProtein = aLine.strip()
  #      ids.append(aProtein)





# http://www.genome.jp/kegg-bin/show_pathway?lif00010+LINJ_35_0990



def from_uni_ids_get_kegg_path_marker(uni_ids):

    # 
    from bioservices import KEGG
    k = KEGG(verbose=False)
    # ORGANISM : hsa
    # KEGG ID : hsa1234
    # PATHWAY : hsa12345    

    print "\n\tKEGG  -  Look for relevant Pathways, using ac. . . \n Protein Marker "

    from bioservices.uniprot import UniProt
    u = UniProt(verbose=False)
    kegg_ids = []
    organism_ids = []
    #filename = raw_input("Enter file name : ")
    af = open("kegg_data\\kegg"+str(uni_ids[0])+"_get_path_marker.txt","wt") # +str(filename)+
    af.write("UniProt ID \tKegg Organism \tProtein Name \tKegg Pathways\n") #str(uni_i77d)+"\t"+str(k_organism)+"\t"+str(get_pname[:-2])+"\t"+str(kegg_paths1)+"\t"+str(paths))
    count = 1
    
    one_wait = 11 # seconds 2 = 27, 4=55, 6=1.8~107, 12=3.19~200, 179 =24:46
    expected_wait = len(uni_ids)*one_wait
    x = (expected_wait)/60
    g = float("{0:.1f}".format(x))
    print "\nTotal ids: ",len(uni_ids)," >>>  Expected job time  : ", g, " minutes. . . "
    for uni_id in uni_ids:
        
        print "\n  Mapping ACC to KEGG ID . . . "
        map1  = u.mapping(fr="ACC+ID", to="KEGG_ID", query=str(uni_id)) # AC
        #print len(map1), map1 # 1 defaultdict(<type 'list'>, {u'P43403': [u'hsa:7535']})
        for genes in map1.values():
            #i = str(i) # [u'lif:LINJ_10_0520'] or [u'lif:LINJ_28_2050', u'lif:LINJ_28_2060']
            #print genes # list
            for i in genes:
                #print  i
                kegg_id = i # # lif:LINJ_10_0520
                #print kegg_id
                kegg_ids.append(kegg_id)
                _id = i[4:] # LINJ_10_0520
                k_organism = i[0:3] # 'lif'
                #print uni_id, " > ", kegg_id #  > lif:LINJ_10_0520
                #kegg_ids.append(kegg_id)
        if k_organism not in organism_ids:
            organism_ids.append(k_organism)

        #print k.get(str(kegg_id))
        get_kegg = k.get(str(kegg_id)) # ENTRY (SEQ...)

        get_pname = k.find(str(k_organism), str(kegg_id))
        get_pname = str(get_pname) # NAMES

        print count , uni_id, " > ",get_pname[:-1] # lif:LINJ_14_0700	putative fatty acid elongase (EC:2.3.1.119)
        
            

        #print "Find pathway from gene . . . "
        #print help(map1)
        #print map1.values()
        #i = str(i)
        # print i

        get_path = k.get_pathway_by_gene(str(_id), str(k_organism)) # PATH
        # if no path is found, will print "No pathway found ?"
        
        #print "! No pathway found ! "
        #print get_path # if empty, will print None
        #paths = str(get_path)
        #dict2 = eval(str1)
        #print dict1==dict2
        #print type(get_path)
            
        
        if get_path.__repr__() == "None":
            #print get_path
            kegg_path = "None"
            af.write(str(uni_id)+"\t"+str(k_organism)+"\t"+str(get_pname[:-1])+"\t"+str(kegg_path)+"\n")

        elif isinstance(get_path, unicode):# isinstance('hello', str)
            #get_path ={}
            pathways = []
            
            kegg_path = get_path[0:8]
            print kegg_path
            kegg_path_name = get_path[10:]
            print kegg_path_name
            kegg_name_fix = kegg_path_name.replace("/","-")
            

            pathways.append(kegg_path+":"+kegg_path_name)
            cf = open("kegg_data\\"+str(uni_id)+"_"+str(kegg_path)+"_"+str(kegg_name_fix)+"_path_mark.png","wb") # kegg_path_name
            res1 = []

            try:
                # A4IAZ4
                print "\nSaving image of the pathway  ",kegg_path, kegg_path_name
                #show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
                res1 =  k.get(str(kegg_path)+"/image") # same as : res =  s.get("hsa05130","image")
                cf.write(res1) # IMAGE
                #for i in pathways:
                print "Showing Pathway Marker in internet page . . . "
                #ef = open(str(kegg_path)+"_"+str(kegg_id)+"_path.png","wb")
                show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})

            except:
                cf.close()

            kegg_paths = " |".join(pathways)
            print "\n> Pathway found: ", kegg_paths
            line = (str(uni_id)+"\t"+str(k_organism)+"\t"+str(get_pname[:-1])+"\t"+str(kegg_paths))
            line.replace("\n","\t")
            af.write(line+"\n")

            #for i in pathways:
            #print "Showing Pathway Marker in internet page . . . "
            #ef = open(str(kegg_path)+"_"+str(kegg_id)+"_path.png","wb")
            #show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
                
            

        else:
            items = []
            pathways = []
            kegg_paths = []
            for i in get_path.values():
                #print i
                items.append(i)
            words = " |".join(items)
            #bf.write(str(words))
            for i in get_path.keys():
                if i.isdigit() == True and len(i) == 5:
                    #pathways.append(str(i)+":"+str(get_path[i]))# = str(i)
                    #print k_organism, i, map1[i]
                    kegg_path = i
                    kegg_path = k_organism + kegg_path
                    kegg_name_fix1 = get_path[i]
                    kegg_name_fix = kegg_name_fix1.replace("/","-")
                    
                    pathways.append(kegg_path+":"+get_path[i])
                    # cf = open("kegg_data\\"+str(uni_id)+"_"+str(kegg_path)+"_"+str(kegg_name_fix)+"_path.png","wb") # 
                    namex = "kegg_data\\"+str(kegg_path)+"_"+str(kegg_name_fix)+"_path_mark.png"
                    if os.path.isfile(namex) == False:
                        cf = open(namex,"wb")
                        res1 = []
                        try:
                            # A4IAZ4
                            print "\nSaving image of the pathway  ",kegg_path, get_path[i]
                            #show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
                            res1 =  k.get(str(kegg_path)+"/image") # same as : res =  s.get("hsa05130","image")
                            cf.write(res1) # IMAGE
                            #for i in pathways:
                            print "Showing Pathway Marker in internet page . . . "
                            #ef = open(str(kegg_path)+"_"+str(kegg_id)+"_path.png","wb")
                            show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
                        except:
                            cf.close()
                        
            #ef.write(str(show_path))
            kegg_paths = " |".join(pathways)
            print "\n> Pathway found: ", kegg_paths
            line = (str(uni_id)+"\t"+str(k_organism)+"\t"+str(get_pname[:-1])+"\t"+str(kegg_paths))
            #line.replace("\n","\t")
            af.write(line+"\n")

            # res0 =  k.get("hsa05130/image") # NF-kappa B signaling pathway - Homo sapiens (human)

            #for i in pathways:
            print "Showing Pathway Marker in internet page . . . "
            #ef = open(str(kegg_path)+"_"+str(kegg_id)+"_path.png","wb")
            show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
            
                

        #print k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
        bf = open("kegg_data\\"+str(uni_id)+"_kegg_entry.txt","w")
        
        bf.write(get_pname+"\n") # NAME
        bf.write(get_kegg) # ENTRY

        count +=1
        

    af.close()
    bf.close()
    print af, bf,"\n", " saved!"
    
    for i in organism_ids:
        df = open("kegg_data\\"+str(i)+"_kegg_pathways.txt","wt")
        #print k.list("pathway", organism= i)
        res_paths = k.list("pathway", organism= i)
        pathways = [x.split()[0] for x in res_paths.strip().split("\n")]
        print len(pathways) ,"pathways found in Kegg organism ", i # as of Nov 2014 -> 94
        df.write(res_paths)
        df.close()
        print df, "saved!"
        
    

    #map2 = k.mapping(fr='KEGG_ID', to='ACC', format='tab', query='hsa:7535') # to="PDB_ID"
    #for i in map2.values():
    #print i#'P43403'
    print "Total ACS : ", count-1
    print "\nBioservices - Kegg - ACCESS ACS - DONE!"
    theInputFiles.append(df.name)
    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    print "(complete!)"


#uni_ids = ["A4HVU4","A4I2W5"] # ["A4I0R1","A4HYI9"] #"A4I5T4","A4I0N4"]#"A4I9M5","A4HVK7","E9AHP1","A4HZH7","A4I8Y5","E9AI05","A4HRR9","A4I5T4","A4I0N4","A4I9M5","A4HVK7","E9AHP1","A4HZH7","A4I8Y5","E9AI05","A4HRR9"]

#from_uni_ids_get_kegg_path_marker(uni_ids)

    

#######################################################################################
    
def from_uni_id_get_kegg_orthologues(uni_id):
    
    import webbrowser
    from bioservices.uniprot import UniProt

    
    # ORGANISM : hsa
    # KEGG ID : hsa1234
    # PATHWAY : hsa12345

    k = KEGG(verbose=False)

    print "\nFrom gene get kegg orthologues . . . "

    from bioservices.uniprot import UniProt
    u = UniProt(verbose=False)

    kegg_ids = []
    print "Mapping protein ac to kegg ID. . . "
    # http://pythonhosted.org/bioservices/convertor_tutorial.html
    map1  = u.mapping(fr="ACC+ID", to="KEGG_ID", query=str(uni_id))# P_ENTREZGENEID, "UNIGENE_ID, EMBL,GeneID
    print map1 # defaultdict(<type 'list'>, {u'P43403': [u'hsa:7535']})

    for i in map1.values():
        i = str(i)
        kegg_id = i[3:-2]
        _id = i[7:-2]
        print _id # LINJ_10_0520
        print i, " = ", kegg_id # [u'lif:LINJ_10_0520'] = lif:LINJ_10_0520
        k_organism = i[3:6] # 'lif'
        print k_organism
        kegg_ids.append(kegg_id)
    print k.get(str(kegg_id))
    get_kegg = k.get(str(kegg_id))

    print "\nOpening page in browser . . ."    
    webbrowser.open("http://www.kegg.jp/ssdb-bin/ssdb_best?org_gene="+str(kegg_id))

    print "BioServices - KEGG access DONE!"
    

#ids = ["A4HUD2","A4HUD4"] 
#for i in ids:
 #   from_uni_id_get_kegg_orthologues(i)



###############################################
 
def from_gene_id_get_kegg_orthologues(gene_id):
    
    import webbrowser
    from bioservices.uniprot import UniProt
    
    letters= ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y","Z"]
    numbers = ["1","2","3","4","5","6","7","8","9","0"]
    
    # ORGANISM : hsa
    # KEGG ID : hsa1234
    # PATHWAY : hsa12345

    k = KEGG(verbose=False)

    print "\nFrom gene get kegg orthologues . . . "

    print "Mapping gene to protein ac . . . "
    u = UniProt()
    res = u.search(str(gene_id), frmt="tab", columns="genes, id")
    for i in res.split("\n"):
        print i
        while len(i.split("\t")) > 1 :
            i = i.split()[-1]
            #print i
            if i[0] in letters and i[1] in numbers and i[-1] in numbers:
                print i
                uni_id = i

    kegg_ids = []

    print "Mapping protein ac to kegg ID. . . "
    # http://pythonhosted.org/bioservices/convertor_tutorial.html
    map1  = u.mapping(fr="ACC+ID", to="KEGG_ID", query=str(uni_id))# P_ENTREZGENEID, "UNIGENE_ID, EMBL,GeneID
    print map1 # defaultdict(<type 'list'>, {u'P43403': [u'hsa:7535']})

    for i in map1.values():
        i = str(i)
        kegg_id = i[3:-2]
        _id = i[7:-2]
        print _id # LINJ_10_0520
        print i, " = ", kegg_id # [u'lif:LINJ_10_0520'] = lif:LINJ_10_0520
        k_organism = i[3:6] # 'lif'
        print k_organism
        kegg_ids.append(kegg_id)
    print k.get(str(kegg_id))
    get_kegg = k.get(str(kegg_id))

    print "\nOpening page in browser . . ."    
    webbrowser.open("http://www.kegg.jp/ssdb-bin/ssdb_best?org_gene="+str(kegg_id))

    print "BioServices - KEGG access DONE!"

#ids = ["LINJ_10_0240","LINJ_10_0520"] 
#for i in ids:
 #   from_gene_id_get_kegg_orthologues(i)

#################################
 
def from_uni_id_get_kegg_spot(uni_id):
    
    from bioservices import KEGG
    k = KEGG(verbose=False)
    # ORGANISM : hsa
    # KEGG ID : hsa1234
    # PATHWAY : hsa12345    

    print "\n\tKEGG  -  Look for relevante Pathways, using ac. . . "
    from bioservices.uniprot import UniProt
    u = UniProt(verbose=False)
    map1  = u.mapping(fr="ACC", to="KEGG_ID", query=str(uni_id))
    print map1 # defaultdict(<type 'list'>, {u'P43403': [u'hsa:7535']})

    kegg_ids = []
    for i in map1.values():
        i = str(i)
        kegg_id = i[3:-2]
        _id = i[7:-2]
        #print _id # lif:LINJ_10_0520
        print i, " = ", kegg_id # [u'lif:LINJ_10_0520'] = lif:LINJ_10_0520
        k_organism = i[3:6] # 'lif'
        #print k_organism
        kegg_ids.append(kegg_id)
        print k.get(str(kegg_id))

    get_path = None
    print k.list("pathway", organism= k_organism)
    res_paths = k.list("pathway", organism= k_organism)
    pathways = [x.split()[0] for x in res_paths.strip().split("\n")]
    print len(pathways) ,"pathways found in Kegg organism ", k_organism # as of Nov 2014 -> 94

  
    print uni_id, "\n",k.find(str(k_organism), str(kegg_id))
    get_pname = k.find(str(k_organism), str(kegg_id))

    

    print "Find pathway from gene . . . "
    print k.get_pathway_by_gene(str(_id), str(k_organism))
    get_path = k.get_pathway_by_gene(str(_id), str(k_organism)) #else "null"

    #bf = open(str(kegg_path)+"_path.txt","w")
    bf = open(str(uni_id)+"_"+str(k_organism)+"_pathways.txt","w")


    if get_path == None:
        print "NO PATHWAY FOUND FOR : ", uni_id
        
        bf.write(res_paths+"\n")
        bf.write(get_pname)
        # break
    else:
        bf.write(res_paths+"\n")
        bf.write(get_pname)
        for i in get_path.keys():
            i = str(i)
            bf.write(str(i)+"\n")
            if i.isdigit() == True and len(i) == 5:
                print i
                kegg_path = i
                kegg_path = k_organism + kegg_path
                print "\n",uni_id, " = ",kegg_id, "\nPathway found: ", kegg_path
            elif len(i) == 8:
                
                kegg_path = i
                print "\n",uni_id, " = ",kegg_id, "\nPathway found: ", kegg_path

            print "Showing Pathway Marker in internet page . . ."
            #print k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
            show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})

            

            "Save image of the pathway : "
            # NF-kappa B signaling pathway - Homo sapiens (human)
            # res0 =  k.get("hsa05130/image")
            res1 =  k.get(str(kegg_path)+"/image") # same as : res =  s.get("hsa05130","image")
            cf = open(str(uni_id)+"_path.png","wb")
            cf.write(res1)
        bf.close()
        cf.close()
        print bf,"saved ! \n",cf, " saved!"
        theInputFiles.append(bf.name)
            

    #map2 = k.mapping(fr='KEGG_ID', to='ACC', format='tab', query='hsa:7535') # to="PDB_ID"
    #for i in map2.values():
    #print i#'P43403'
    print "\nBioservices - Kegg - ACCESS AC - DONE!"
    #print f1, "\n",f2,"\n", f3,"\n", f4 ,"\n", " are saved!"

   
    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    print "(complete!)"


#ids_kegg = ["A4HSF7","A4HZJ3","Q7K8Z6"]
#ids_kegg = ["A4HXU4"]
#ids_kegg = ["A4HXL7","A4IBW5"]#,"A4I599"]
#for i in ids_kegg:
#    from_uni_id_get_kegg(i)


######################################################################################3

def from_gene_id_get_kegg(gene_id):
    
    letters= ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y","Z"]
    numbers = ["1","2","3","4","5","6","7","8","9","0"]
    
    # ORGANISM : hsa
    # KEGG ID : hsa1234
    # PATHWAY : hsa12345

    k = KEGG(verbose=False)

    print "KEGG  -  Look for relevante Pathways, using gene . . . "

    print "Mapping gene to protein ac . . . "
    u = UniProt()
    res = u.search(str(gene_id), frmt="tab", columns="genes, id")
    for i in res.split("\n"):
        print i
        while len(i.split("\t")) > 1 :
            i = i.split()[1]
            #print i
            if i[0] in letters and i[1] in numbers and i[-1] in numbers:
                print i
                uni_id = i

    kegg_ids = []

    print "Mapping protein ac to kegg ID. . . "
    # http://pythonhosted.org/bioservices/convertor_tutorial.html
    map1  = u.mapping(fr="ACC+ID", to="KEGG_ID", query=str(uni_id))# P_ENTREZGENEID, "UNIGENE_ID, EMBL,GeneID
    print map1 # defaultdict(<type 'list'>, {u'P43403': [u'hsa:7535']})

    for i in map1.values():
        i = str(i)
        kegg_id = i[3:-2]
        _id = i[7:-2]
        print _id # LINJ_10_0520
        print i, " = ", kegg_id # [u'lif:LINJ_10_0520'] = lif:LINJ_10_0520
        k_organism = i[3:6] # 'lif'
        print k_organism
        kegg_ids.append(kegg_id)
    print k.get(str(kegg_id))
    get_kegg = k.get(str(kegg_id))


    print "Find  protein . . . "
    #pname = raw_input("Please enter protein name : ")
    print k.find(str(k_organism), uni_id)
    get_pname = k.find(str(k_organism), uni_id)

    print "Find  pathways for this organism . . . "
    print k.list("pathway", organism= k_organism)
    res_paths = k.list("pathway", organism= k_organism)
    pathways = [x.split()[0] for x in res_paths.strip().split("\n")]
    print len(pathways) ,"pathways found in Kegg organism ", k_organism # as of Nov 2014 -> 94
    
   
    print uni_id, "\n",k.find(str(k_organism), str(kegg_id)) # lif:LINJ_14_0700	putative fatty acid elongase (EC:2.3.1.119)

    get_pname = k.find(str(k_organism), str(kegg_id))

    print "Find pathway from gene . . . "
    get_path = k.get_pathway_by_gene(str(_id), str(k_organism))
    
    if get_path.__repr__() == "None":
        print get_path
        kegg_path = "None"
    else:
        for i in get_path.keys():
            print i
            i = str(i)
            if i.isdigit() == True and len(i) == 5:
                print i
                kegg_path = i
                kegg_path = k_organism + kegg_path
                print "\n",uni_id, " = ",kegg_id, "\nPathway found: ", kegg_path

    print "Showing Pathway Marker in internet page . . ."
    #print k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
    
    bf = open(str(gene_id)+"_path.txt","w")
    bf.write(res_paths+"\n")
    bf.write(get_pname+"\n")
    bf.write(get_kegg)

    if get_path.__repr__() != "None":
        for i in get_path:
            bf.write(str(i))
    
    print "Save image of the pathway . . . "
    # NF-kappa B signaling pathway - Homo sapiens (human)
    # res0 =  k.get("hsa05130/image")
    if kegg_path != "None":
        show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
        res1 =  k.get(str(kegg_path)+"/image") # same as : res =  s.get("hsa05130","image")
        cf = open(str(uni_id)+"_path.png","wb")
        cf.write(res1)
    cf.close()
    print cf, " saved!"

    bf.close()
    print bf,"\n", " saved!"

    #map2 = k.mapping(fr='KEGG_ID', to='ACC', format='tab', query='hsa:7535')
    #for i in map2.values():
    #print i#'P43403'
    print "\nBioservices - Kegg - ACCESS Gene - DONE!"
    #print f1, "\n",f2,"\n", f3,"\n", f4 ,"\n", " are saved!"
    theInputFiles.append(bf.name)
    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    print "(complete!)"

#from_gene_id_get_kegg("LINJ_14_0700")

    
# http://www.genome.jp/dbget-bin/www_bget?lif:LINJ_14_0700
# http://www.ncbi.nlm.nih.gov/gene/?term=LINJ_14_0700

# http://pythonhosted.org/bioservices/_modules/bioservices/uniprot.html
# http://central.biomart.org/martsearch/#!/?q=
#http://biodbnet.abcc.ncifcrf.gov/db/db2db.php#biodb

# http://www.cureffi.org/2013/03/15/mapping-proteomics-data-to-uniprot-refseq-and-gene-symbol/

# http://www.kegg.jp/dbget-bin/www_bget?lif:LINJ_10_0490+lif:LINJ_10_0500+lif:LINJ_10_0510+lif:LINJ_10_0520+lif:LINJ_10_0530+lif:LINJ_28_0600+lif:LINJ_28_0610+lif:LINJ_31_2040
# http://www.kegg.jp/dbget-bin/www_bget?lma:LMJF_10_0465
# https://pythonhosted.org/bioservices/protein.html
#######################################################################################



#######################################################################################
# NCBI
# ((rod) AND "eukaryotes"[porgn:__txid2759]) AND "kinetoplastids"[porgn:__txid5653] 
######################################

# Conseved Domains NCBI
# http://www.ncbi.nlm.nih.gov/Structure/cdd/wrpsb.cgi

##############################################
#####################################################################################################33


def from_uni_ids_get_kegg_path(uni_ids):
        
        # FROM UNI ID CONVERT TO KEGG ID AND SHOWS in PATH IN WEB PAGE (red colour)

   import datetime, time
   time_begin = datetime.datetime.fromtimestamp(time.time())
     
   
   k = KEGG(verbose=False)
   
   #k.lookfor_organism("leish")
   print k.lookfor_organism(raw_input("\nLook for organism kegg name : "))# "droso" 

   #print "KEGG Tax ID: ", k.code2Tnumber("hsa") # T01001
   #print k.get("hsa:7535")


   k.organism = raw_input("\nEnter organism kegg name (3 letters): ")#"lif"
   print "Converting KEGG Tax ID to Tnumber ", k.code2Tnumber(k.organism) # T01001

   #if k.organism == True:
   print k.list("pathway", organism= k.organism)
   res = k.list("pathway", organism= k.organism)
   pathways = [x.split()[0] for x in res.strip().split("\n")]
   print len(pathways) ,"pathways found in Kegg organism ", k.organism # as of Nov 2014 -> 94
    
    #print s.get("lif:LINJ_10_0520")
    # INFO:root:REST.bioservices.Kegg request begins
    # INFO:root:--Fetching url=http://rest.kegg.jp/get/lif:LINJ_08_0960
   
   from bioservices.uniprot import UniProt
   u = UniProt(verbose=False)
   # CONVERT TO KEGG ID
   # u.mapping(fr="ACC+ID", to="KEGG_ID", query='P43403')
   # {'P43403': ['hsa:7535']}
   kegg_ids = []
   for i in uni_ids:
           res_id = u.mapping("ACC", "KEGG_ID", query=str(i))
           kegg_id = res_id.values()
           print ', '.join("%s = %r" % (key,val) for (key,val) in res_id.iteritems())
           print kegg_id           
           #kegg_ids = kegg_ids.append(kegg_id)
           for i in res_id.values():
                   i = str(i)
                   print i
                   kegg_ids.append(i[3:-2])

        # GET ENTRY - hsa:7535
   for i in kegg_ids:
       print i
       print k.get(str(k.organism)+":"+str(i))
   #print k.get("hsa:7535")
   
   
   """
   ENTRY       7535              CDS       T01001
   NAME        ZAP70, SRK, STCD, STD, TZK, ZAP-70
   DEFINITION  zeta-chain (TCR) associated protein kinase 70kDa (EC:2.7.10.2)
   ORTHOLOGY   K07360  tyrosine-protein kinase ZAP-70 [EC:2.7.10.2]
   ORGANISM    hsa  Homo sapiens (human)
   PATHWAY     hsa04014  Ras signaling pathway
               hsa04064  NF-kappa B signaling pathway
               hsa04650  Natural killer cell mediated cytotoxicity
               hsa04660  T cell receptor signaling pathway
               hsa05340  Primary immunodeficiency
   DISEASE     H00093  Combined immunodeficiencies (CIDs)
   BRITE       KEGG Orthology (KO) [BR:hsa00001]
   ...
   DBLINKS : UniProt: P43403
   """
   
   print "Find gene from name: (organism, protein_name)"
   pname = raw_input("\nEnter protein name : ")
   print k.find(k.organism, pname)

   # GET PATHWAY BY GENE OR KEGGID + k.organism
   # Use the correspondence number to search kegg pathways 
   # Get pathways that contain this KEGG ID:
   """for i in kegg_ids:
           #print k.get_pathway_by_gene("7535", "hsa") , "\n"
           print k.get_pathway_by_gene("LINJ_08_0960", "lif")
           res= k.get_pathway_by_gene("LINJ_08_0960", "lif")
           print ', '.join("%s = %r" % (key,val) for (key,val) in res_id.iteritems())
           res_paths=[]
           for i in res.values():
                   i = str(i)
                   print i
                   res_paths.append(i[3:-2])


           
           
           #print k.get_pathway_by_gene(str(i), "lif")
           print k.find("pathway", str(i))
           #lif05140
           # LINJ_30_0420
           print "see" 
      """     
   """{'hsa04660': 'T cell receptor signaling pathway', 'hsa04650': 'Natural killer cell mediated cytotoxicity',
   'hsa04064': 'NF-kappa B signaling pathway', 'hsa04014': ' Ras signaling pathway',
   'hsa05340': 'Primary immunodeficiency'}"""
   
   #####################  Show in INTERNET page  ##############
   
   #print k.show_pathway("hsa04064", keggid={"7535": "red"})
   kegg_id1 = "LINJ_08_0950"
   print k.show_pathway("lif05140", keggid={str(kegg_id1): "red"}) #CBP in Prostate Cancer
   
   # http://www.genome.jp/kegg-bin/show_pathway?map=hsa05215&show_description=show
   # http://www.genome.jp/kegg-bin/show_pathway?map=lif05140&show_description=show
   
   # http://www.genome.jp/kegg-bin/show_pathway?org_name=lif&mapno=00061&mapscale=&show_description=show
   # http://www.genome.jp/kegg-bin/show_pathway?map=lif00061
   
   # http://www.genome.jp/dbget-bin/www_bget?hsa:1387+hsa:2033
   # http://www.genome.jp/dbget-bin/www_bget?lif:LINJ_24_1560
   
   ### Save image of the pathway :
   # NF-kappa B signaling pathway - Homo sapiens (human)
   # res0 =  k.get("hsa05130/image")
   ###res1 =  k.get("lif05140/image")
   # same as : res =  s.get("hsa05130","image")
   k.show_pathway(str(k_organism)+str(keggg_id), keggid={str(kegg_id): "red"})
   f = open(str(kegg_id1)+".png", "wb")
   f.write(res1)
   f.close()
   #----------------------
   time_end = datetime.datetime.fromtimestamp(time.time())
   print("Time elapsed: ", str(time_end - time_begin))

#uni_ids= ["A4HUG0"," "]
#from_uni_ids_get_kegg_path("A4HUG0")





#################################################33

def from_gene_ids_get_kegg_path(gene_ids):

    s = KEGG()
    #s.lookfor_organism("droso")
    #['T01014 lma Leishmania major Eukaryotes;Protists;Euglenozoa;Kinetoplasts', 'T01112 lif Leishmania infantum Eukaryotes;Protists;Euglenozoa;Kinetoplasts', 'T02289 ldo Leishmania donovani Eukaryotes;Protists;Euglenozoa;Kinetoplasts', 'T02288 lmi Leishmania mexicana Eukaryotes;Protists;Euglenozoa;Kinetoplasts', 'T01113 lbz Leishmania braziliensis Eukaryotes;Protists;Euglenozoa;Kinetoplasts']
    #print s.find("pathway", "B+cell")
    #print s.lookfor_pathway("B cell")
    s.organism = raw_input("Enter organism kegg name : ")#"lif"
    
    if s.organism == True:
            print s.list("pathway", organism= s.organism)
            res = s.list("pathway", organism= s.organism)
            pathways = [x.split()[0] for x in res.strip().split("\n")]
            print len(pathways) ,"pathways" # as of Nov 2014 -> 94
    
    #print s.get("lif:LINJ_10_0520")
    # INFO:root:REST.bioservices.Kegg request begins
    # INFO:root:--Fetching url=http://rest.kegg.jp/get/lif:LINJ_08_0960
    """
    ENTRY       LINJ_08_0960      CDS       T01112
    
    DEFINITION  cathepsin L-like protease
    ORTHOLOGY   K13537  cysteine peptidase B [EC:3.4.22.-]
    ORGANISM    lif  Leishmania infantum
    PATHWAY     lif05140  Leishmaniasis
    BRITE       KEGG Orthology (KO) [BR:lif00001]
                 Human Diseases
                 ...
    POSITION    8:complement(410134..411465)
    MOTIF       Pfam: Peptidase_C1 DUF3586 Inhibitor_I29 Peptidase_C1_2
    DBLINKS     NCBI-GI: 339896953
                NCBI-GeneID: 10966317
                UniProt: Q7K8Z6
    AASEQ       443
                MATSRAALCAVAVVCVVLAAACAPARAIYVGTPAAALFEEFKRTYRRAYGTLAEEQQRLA
                NFERNLELMREHQARNPHARFGITKFFDLSEAEFAARYLNGAAYFAAAKQHAGQHYRKAR
                ADLSAVPDAVDWREKGAVTPVKNQGACGSCWAFSAVGNIESQWARAGHGLVSLSEQQLVS
                CDDKDNGCNGGLMLQAFEWLLRHMYGIVFTEKSYPYTSGNGDVAECLNSSKLVPGAQIDG
                YVMIPSNETVMAAWLAENGPIAIAVDASSFMSYQSGVLTSCAGDALNHGVLLVGYNKTGG
                VPYWVIKNSWGEDWGEKGYVRVVMGLNACLLSEYPVSAHVPRSLTPGPGTESEERAPKRV
                TVEQMMCTDMYCREGCKKSLLTANVCYKNGGGGSSMTKCGPQKVLMCSYSNPHCFGPGLC
                LETPDGKCAPYFLGSIMNTCQYT"""
    #print s.get("uniprot:P43403")
    print s.get("path:"+str(s.organism)+raw_input("Enter Path ID (5 digits) :"))
    # INFO:root:REST.bioservices.Kegg request begins
    # INFO:root:--Fetching url=http://rest.kegg.jp/get/path:lif01100
    """
    ENTRY       lif01100                    Pathway
    NAME        Metabolic pathways - Leishmania infantum
    PATHWAY_MAP lif01100  Metabolic pathways
    ORGANISM    Leishmania infantum [GN:lif]
    KO_PATHWAY  ko01100
    ///"""
    #print s.find("lif", "PKC")                          NAME
    #print s.find("lif", "ko01100")                      CODE
    for gene in gene_ids:
            print s.get_pathway_by_gene(str(gene), s.organism)#  GENE
            #print s.get_pathway_by_gene(gene_id, s.organism)#  GENE

    print "\nBioservices - Kegg - ACCESS Gene - DONE!"        
    #print f1, "\n",f2,"\n", f3,"\n", f4 ,"\n", " are saved!"

    

#from_gene_id_get_kegg("LINJ_14_0700")
    

##################################3

def from_kegg_ids_get_kegg_path(kegg_ids):

    # FROM UNI ID CONVERT TO KEGG ID AND SHOWS in PATH IN WEB PAGE (red colour)

    import datetime, time
    time_begin = datetime.datetime.fromtimestamp(time.time())

    k = KEGG(verbose=False)

    # CONVERT TO KEGG ID
    from bioservices.uniprot import UniProt
    u = UniProt(verbose=False)

    # u.mapping(fr="ACC+ID", to="KEGG_ID", query='P43403')
    # {'P43403': ['hsa:7535']}

    for i in uni_ids:
        if len(i)<6:
                pass
        res_id = u.mapping("ACC", "KEGG_ID", query=str(i))
        kegg_id = res_id.values()
        print i , " > ", kegg_id
        #print ', '.join("%s = %r" % (key,val) for (key,val) in res_id.iteritems())

        kegg_ids = []
        for i in res_id.values():
            i = str(i)
            print i
            kegg_ids.append(i[3:-2])

    # GET ENTRY - hsa:7535
    for i in kegg_ids:
        print i
        #print k.get(str(k.organism)+":"+str(i))

    
#uni_ids= ["A4HUG0"," "]
#from_uni_ids_get_kegg_path(uni_ids)

def from_gene_uni_ids_get_kegg_image():
    print "In Progress . . ."

######################################################################################
####################33             Look for pathways (by name)          ###################

def from_name_get_kegg_paths_genes():
        #  https://pythonhosted.org/bioservices/kegg_tutorial.html

        # FROM UNI ID CONVERT TO KEGG ID AND SHOWS in PATH IN WEB PAGE (red colour)


        import datetime, time
        time_begin = datetime.datetime.fromtimestamp(time.time())

        from bioservices import KEGG
        k = KEGG(verbose=False)

        #k.lookfor_organism("leish")
        print "\nLook for pathways (by name)"
        print k.lookfor_organism(raw_input("\nLook for organism kegg name : "))# "droso"

        #print "KEGG Tax ID: ", k.code2Tnumber("hsa") # T01001
        #print k.get("hsa:7535")

        k.organism = raw_input("\nEnter organism kegg name (3 letters): ")
        print "Converting KEGG Tax ID to Tnumber ", k.code2Tnumber(k.organism) # T01001
        #if k.organism == True:

        print k.list("pathway", organism= k.organism)

        res = k.list("pathway", organism= k.organism)
        pathways = [x.split()[0] for x in res.strip().split("\n")]
        print len(pathways) ,"pathways found in Kegg organism ", k.organism # as of Nov 2014 -> 94

        f1 = open(str(k.organism)+"_pathways.txt","wt")
        for i in res.strip().split("\n"):
                f1.write(i+"\n")
        f1.write("Total pathways : "+str(len(pathways)))
        f1.close()


        print "\nFind gene from protein name : (organism, protein_name)"
        name = raw_input("\nEnter protein name : ")
        f2 = open(str(k.organism)+"_"+name+"_genes.txt","wt")
        for i in k.find(k.organism, name).splitlines():
                print i
                f2.write(i+"\n")

        f3 = open(str(k.organism)+"_"+name+"_path.txt","wt")
        f4 = open(name+".png", "wb")

        # USE GENE
        print k.get(str(k.organism)+":"+"LINJ_05_0040")
        for i in k.get_pathway_by_gene("LINJ_08_0960", k.organism):
                if i[0:3] == k.organism:
                        print i
                        print k.get(i)
                        f3.write(k.get(i))
                        # Show in INTERNET page
                        print k.show_pathway(i)
                        res1 =  k.get(str(i)+"/image")
                        f4.write(res1)
        
        f2.close()
        f3.close()        
        f4.close()

        # Show in INTERNET page
        # print k.show_pathway("hsa04064", keggid={"7535": "red"})
        # print k.show_pathway("lif05140", keggid={"LINJ_08_0950": "red"}) #CBP in Prostate Cancer

        # Save image of the pathway :
        # NF-kappa B signaling pathway - Homo sapiens (human)
        # res0 =  k.get("hsa05130/image")

        
        # same as : res =  s.get("hsa05130","image")

        print "\nBioservices - Kegg - ACCESS Name - DONE!"        
        print f1, "\n",f2,"\n", f3,"\n", f4 ,"\n", " are saved!"

       
  

#from_name_get_kegg_paths_genes()


################################################################################
##### EBI                                                                 ######
################################################################################

# http://www.ebi.ac.uk/interpro/download.html




################################################################################
##### PDB                                                                 ######
################################################################################
#from bioservices import PDB

def from_pdb_id_get_pdb(pdb_id):

    from bioservices import PDB
    
    s = PDB()
    # for i in dir(s): print i
    """['CACHE_NAME', 'CACHING', 'TIMEOUT', '__class__', '__copy__', '__deepcopy__', '__delattr__', '__dict__', '__doc__', '__format__',
        '__getattribute__', '__hash__', '__init__', '__module__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__',
        '__sizeof__', '__str__', '__subclasshook__', '__wakref__', '_apply', '_create_cache_session', '_create_session', '_debugLevel',
        '_easyXMLConversion', '_get_all_urls', '_get_async', '_get_caching', '_get_easyXMLConversion', '_get_level', '_get_session',
        '_get_timeout', '_get_url', '_interpret_returned_request', '_process_get_request', '_service', '_session', '_set_caching',
        '_set_easyXMLConversion', '_set_level', '_set_timeout', '_set_url', '_url', 'clear_cache', 'content_types', 'critical', 'debug',
        'debugLevel', 'debug_message', 'delete_cache', 'devtools', 'easyXML', 'easyXMLConversion', 'error', 'getUserAgent', 'get_async',
        'get_current_ids', 'get_file', 'get_headers', 'get_ligands', 'get_one', 'get_sync', 'get_xml_query', 'http_delete', 'http_get',
        'http_post', 'http_put', 'info', 'last_response', 'level', 'logging', 'name', 'on_web', 'post_one', 'pubmed', 'requests_per_sec',
        'response_codes', 'save_str_to_image', 'search', 'session', 'settings', 'url', 'warning']
        """
    res1 = s.get_file(pdb_id, "pdb")
    #s.get_current_ids() error

    res2 = s.get_ligands(str(pdb_id))
    print "\n    Bioservices Webservices\nLigand info for ",pdb_id,"\n",res2
    
    #import tempfile
    #fh = tempfile.NamedTemporaryFile()
    #fh.write(res0)
    # manipulate the PDB file with your favorite tool
    # close the file ONLY when finished (this is temporary file)
    #fh.close()

    f1 = open("output/"+str(pdb_id)+".pdb", "w")
    f1.write(res1)
    f1.close()
    print "\nBioservices - PDB - ACCESS PDB STRUCTURE - DONE!"
    # print "Proteins found for this GO term:", count
    print f1 , "is saved!"
    #theInputFiles.append(f1.name)
    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    print "(complete!)"


#from_pdb_id_get_pdb("1FBV")



################################################################################
##### WikiPathways                                                        ######
################################################################################

#import access_quickgo_with_biomart
def access_wikipatways_by_prot_ids(theProtIDs):
        from bioservices import wikipathway
        s= Wikipathway(verbose=False)
        s.getPathwayInfo("WP2320")
        #s = Wikipathway()
        s.organism  # default organism
        #'Homo sapiens'
        s.findPathwaysByText('MTOR')
        s.getPathway('WP1471')
        s.getPathwaysByOntologyTerm('DOID:344')
        s.findPathwaysByXref('P45985')
        for i in theProtIDs:
                s.findPathwaysByXref(str(i))
        print "done!"



#access_wikipatways_by_prot_ids(["O15673","Q94593","O00912","O00913"])
                
        
################################################################################
#####  AUXILIARY FUNCTIONS       - SUMMARY ANALYSIS                                              ######
################################################################################

#def access_kegg_by_gene_name(theGeneList1, theKeegPathways):
        
#  Task 05.02 - QuickGO
#import access_quickgo_with_biomart

#  Task 05.03 - KEGG
#import access_kegg_with_biomart

#############################################################################
#############################################################################
############################################################################


"""Entry	Entry name	Status	Protein names	Gene names	Gene ontology 

IDs	Protein existence	Length	# Proteins	# Unique Peptides	

# Peptides	# PSMs	Area	Medium/Light	Medium/Light Count	

Medium/Light Variability [%]	MW [kDa]	calc. pI
A4IE50	A4IE50_LEIIN	unreviewed	Uncharacterized protein	LINJ_36_4440	

	Predicted	1123	1	11	11	18	3.368E8	

100.000	2	0.0	122.7	5.52
"""

# This function opens a file A and file B to match keys,
# The output is a tab grid

# SEARCH by
# 1) tax id >>> Bioservices
# 2) Upload Files:
#       2.1)sequence
#       2.2)uniprot id
#       2.3)domain
#       2.4)
# 



#######################################################################
##################################################3

from bioservices import *

def get_bmq_ipr2go():

    s = BioMart()
    ret = s.registry() # to get information about existing services aka databases
    #dir(s)
    s = BioMart(verbose=False)
    #print "shadow"
    #xmlq = s.biomartQuery.get_xml()
    #res = s.query(xmlq)
    #s.new_query()
    s.lookfor("interpro")                           #Candidate:
    s.datasets("prod-intermart_1")                  # MART name:)

    ### s.add_dataset_to_xml("protein")
    ### s.add_attribute_to_xml("protein_accession")
    s.add_dataset_to_xml("entry")
    #s.add_attribute_to_xml("protein_ac")
    s.add_attribute_to_xml("entry_id")#entry_name
    s.add_attribute_to_xml("entry_name")
    s.add_attribute_to_xml("entry_type")
    s.add_attribute_to_xml("go_id")
    s.add_attribute_to_xml("go_term_name")
    s.add_attribute_to_xml("go_root_term")
    #s.add_attribute_to_xml("database_name")

    #s.add_filter_to_xml("protein_tax_id_filter", str(tax_id) )
    xml_query = s.get_xml()
    res2 = s.query(xml_query)
    res2 = res2.split("\n") # TOTAL = 42128, 62706(6out)
    header22 ="""entry_id\tentry_name\tentry_type\tgo_id\tgo_term_name\tgo_root_name\n"""
    print header22
    for i in res2[0:10]: print i
    bFile = open("output/_bmq_ipr_go.txt", "wt")
    bFile.write(str(header22))
    for aLine in res2:
        bFile.writelines(aLine+"\n")
    bFile.close()
    print  """Biomart Query | Hits : """, str(len(res2))
    print bFile, " saved!"
    
    
def FileCheck(kf):
    try:
      open(kf, "r")
      return 1
    except IOError:
      print "Error: File %s not appear to exist.\nLet 's create <_bmq_ipr_go.txt>" % kf
      get_bmq_ipr2go()
      return 0

#FileCheck("output/_bmq_ipr_go.txt")
#get_bmq_ipr2go()

#result = FileCheck("testfile")
#print result
    



#####################################################################

def get_db_interpro2go():

    webbrowser.open("ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go")
    webbrowser.open("ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt")


def get_webinfo():
    import urllib
    f = urllib.urlopen("http://api.bitcoincharts.com/v1/trades.csv?symbol=mtgoxUSD")
    print f.read()

#get_webinfo()
    
def get_text_from_web(url):

    import requests
    #symbol = "mtgoxUSD"
    #url = 'http://api.bitcoincharts.com/v1/trades.csv?symbol={}'.format(symbol)
    data = requests.get(url)

    # dump resulting text to file
    #with open("trades_{}.csv".format(symbol), "w") as out_f:
    fname = url.split("/")[-1]
    fname = fname[:-4]
    with open(str(fname)+".txt", "w") as out_f:
        out_f.write(data.text)

    #return out_f
    
#get_text_from_web("ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go")


#################################################################################
def update_family_tree():
   
   import urllib2
   url2 = "ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt"
   response = urllib2.urlopen(url2)
   webContent = response.read()
   print webContent[0:300]
   fname = url2.split("/")[-1]
   fname = fname[:-4]
   with open(str(fname)+".txt","w") as out_f:
      out_f.write(webContent)

#update_family_tree()
      
############3##
def copy_to_file():
   # open("out1.txt", "w").writelines([l for l in open("in.txt").readlines() if "tests/file/myword" in l])
   with open("in.txt") as f:
      with open("out.txt", "w") as f1:
         for line in f:
            if "ROW" in line:
               f1.write(line)
               

#########################3
               
from collections import *
def child2family_aux(sourceFile, filterFile, newFile):

    

    source_file = open(sourceFile, "rt")
    filter_file = open(filterFile, "rt")
    filter_list = []
    ipr2family = defaultlist(list)
    for i in filter_file:
                filter_list.append(i.strip())
    new_file = open( newFile , 'w', 0)
    print >> new_file, 'Domain\tDomain Family\tLevel'
    temp_family = ''
    for line in source_file:
        if line[0:2] != '--':
            temp_family = line.strip()
        elif line[0:2] == '--':
            level = 0
            new_child = line.strip()
            while True:
                new_child = new_child[2:]
                level += 1
                if new_child[0:2] != '--':
                    break
            if new_child.split('::')[0] in filter_list:
                print >> new_file, '%s\t%s\t%s' %(new_child, temp_family, level)
                ipr2family[new_child.split('::')[0]]= temp_family

    return ipr2family

sourcefile = "familytree_ipr2go.txt"
filterFile = "2pep_246_down_filter_246_tax5671_bs_uniprot.txt"
newFile = "tree_ipr_test.txt"

#####################






###############################################################
time_end = datetime.datetime.fromtimestamp(time.time())
print("Time elapsed: ", str(time_end - time_begin))

