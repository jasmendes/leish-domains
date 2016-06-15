################################################################################
##### Loading libraries                                                    #####
################################################################################
import sys, os
#print "FILE PATH:", sys.argv[0]
#print "FILE NAME:", os.path.basename(sys.argv[0])
#print sys.path
sys.path.append("C:\Python273\Lib\site-packages")
sys.path.append("C:\Python276\Lib\site-packages")


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
from bioservices import QuickGO
from math import pi
import urllib

Entrez.email = 'fc34880@alunos.fc.ul.pt'
#import os
#os.environ["http_proxy"]="http://username:password@proxy.alunos.di.fc.ul.pt:3128"
#os.environ["Path"]+=";blastpath"

global theUniIDs


global theInputFiles



import datetime, time

import collections
from collections import Counter

from bioservices_functions_tk import *

def write_ids_to_file(theUniIDs,theFilename):
        
        
        af = open(theFilename, "w")
        count = 1
        from bioservices.uniprot import UniProt
        u = UniProt(verbose=False)

        for uni_id in theUniIDs:
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
                                #kegg_ids.append(kegg_id)
                                _id = i[4:] # LINJ_10_0520
                                k_organism = i[0:3] # 'lif'
                                #print uni_id, " > ", kegg_id #  > lif:LINJ_10_0520
                                #kegg_ids.append(kegg_id)
                                print _id
                                af.write(str(count)+"\t"+_id+"\t"+uni_id+"\n")
                count += 1
        af.close()
        

def from_uni_id_get_bs_uniprot(uni_ids,theFilename):
    
        ### SEARCH BY UNI_ID

        from bioservices.uniprot import UniProt
        u = UniProt(verbose=False)
        """
        Database: Bioservices - UniProt
        Correct values are ['citation', 'clusters', 'comments', 'database', 'domains', 'domain', 'ec', 'id', 'entry name',
        'existence', 'families', 'features', 'genes', 'go', 'go-id', 'interpro', 'interactor', 'keywords', 'keyword-id',
        'last-modified', 'length', 'organism', 'organism-id', 'pathway', 'protein names', 'reviewed', 'score', 'sequence', '3d',
        'subcellular locations', 'taxonomy', 'tools', 'version', 'virus hosts']"""
    
        #s = ",".join(uni_ids)
        count = 1
        aFile = open(theFilename, "wt")
        #aFile.write("id\t  protein names\t interpro\t go-id\t go\t genes\t pathway\t subcellular\t locations\t existence\t organism\t organism-id\n")
        for i in uni_ids:
            data = u.search(i, frmt="tab", limit=1000000000, columns="entry name,id")#" organism-id, entry name, id, protein names, interpro,  genes, go-id, go,  pathway, subcellular locations, existence, organism") # limit seacrh (2march2015) =153
            #print help(data)
            res1 = data.split("\n")
            if count == 0:
                for i in res1[:-1]:
                    #print i
                    aFile.write(str(1)+"\t"+i+"\n")
                    count += 1
            else:
                for i in res1[1:-1]:
                    #print i
                    aFile.write(str(1)+"\t"+i+"\n")
                    count += 1

        #for i in res1:
            #print i
            
        print "Proteins found : ", count-1
        aFile.close()
        print aFile, " saved!"
        #theInputFiles.append(aFile.name)
        
        print "\n\tBioservices - UniProt Access ACS - DONE!\n"

#theUniIDs = ["A4I5T4","A4I0N4","A4I9M5"]#,"A4HVK7","E9AHP1","A4HZH7","A4I8Y5","E9AI05","A4HRR9","A4I5T4","A4I0N4","A4I9M5","A4HVK7","E9AHP1","A4HZH7","A4I8Y5","E9AI05","A4HRR9"]

#theFilename = "1_name_ac.txt"
#from_uni_id_get_bs_uniprot(theUniIDs, theFilename)
        




################################################################################        
##### Protein identifiers and sequences                                   ######
################################################################################
# Reads UniProt accession numbers from a file
# Returns a dictionary with groups of tuples of proteins (name, accession)
# The file should contain one protein per line in the following format:
# <protein_group> <protein_name> <protein_accession_number>
def load_proteins_accessions(theFilename):
    aProteins = {}
    aFile = open(theFilename, "rt")
    for aLine in aFile.readlines():
        aRow = aLine.strip().split("\t")
        if len(aRow)==3:
            if aRow[0] in aProteins:
                aProteins[str(aRow[0])].append({'name':aRow[1], 'uniprot':aRow[2]})
            else:
                aProteins[str(aRow[0])] = [{'name':aRow[1], 'uniprot':aRow[2]}]
        elif len(aRow) == 1:
                aProteins[str(aRow[0])].append({'name':aRow[0], 'uniprot':aRow[0]})
        else:
                
                aProteins[str(aRow[0])].append({'name':aRow[1], 'uniprot':aRow[2]})
                        
    aFile.close()
    return aProteins

# Retrieves a protein sequence from UniProt
def get_proteins_sequences(theProteins):
    for aProtein in theProteins:
        aRecord = get_record(aProtein)
        print aRecord
        print ">sp|"+str(aRecord.accessions[0])+"|"+str(aRecord.entry_name)+" "+str(aRecord.features[0][3])
        print aRecord.sequence

# Retrieve and save all proteins from theProteins in a FASTA file on theFilename
def save_proteins_sequences(theProteins, theFilename):
    aFile = open(theFilename, "wt")
    for aProtein in theProteins:
        aRecord = get_record(aProtein['uniprot'])
        
        if len(aRecord.features) == 0:
                aFile.write(">sp|"+str(aRecord.accessions[0])+"|"+str(aRecord.entry_name)+" "+str(aRecord.entry_name)+"\n")
                aFile.write(str(aRecord.sequence)+"\n")
        else:
                aFile.write(">sp|"+str(aRecord.accessions[0])+"|"+str(aRecord.entry_name)+" "+str(aRecord.features[0][3])+"\n")
                aFile.write(str(aRecord.sequence)+"\n")
                
    aFile.close()

# Reads the aminoacids sequences of proteins from a file
# Returns a dictionary with the sequences indexed by the protein accession number
def load_proteins_sequences(theFilename):
    aFile = open(theFilename, "rt")
    aSequences = {}
    for aSequence in SeqIO.parse(aFile,'fasta'):
        aSequences[str(aSequence.id.split("|")[1])] = str(aSequence.seq)
    aFile.close()
    return aSequences
    
################################################################################        
##### UniProt                                                             ######
################################################################################
# Retrieves a protein record from UniProt
def get_record(theAcessionNumber):
    handle = ExPASy.get_sprot_raw(theAcessionNumber)
    try:
        return SwissProt.read(handle)
    except ValueException:
        print "WARNING: Accession %s not found" % theAcessionNumber
        
################################################################################        
##### Homologues search                                                   ######
################################################################################
# Extract the identifier of each homologue protein (e.g.: CCR1_MOUSE)
def extract_homologues_proteins(theFilename):
    aFile = open(theFilename, "rt")
    aHomologuesProteins = []
    for aSequence in SeqIO.parse(aFile,'fasta'):
        aHomologuesProteins.append(str(aSequence.id.split("\t")[4]))
    aFile.close()
    return aHomologuesProteins
    
def save_homologues(theSequences, theNumHomologues):
    for aSequence in theSequences:
        aFile = open("output/"+str(aSequence)+"_homologues.fasta", "wt")
        aResultHandle = NCBIWWW.qblast("blastp", "swissprot", theSequences[aSequence], hitlist_size=theNumHomologues)
        aBlastRecords = NCBIXML.parse(aResultHandle)
        for aRecord in aBlastRecords:
            for aAlignment in aRecord.alignments:
                aRecord = get_record(aAlignment.accession)
                aFile.write(">"+aAlignment.accession+" mol:protein length:"+str(aAlignment.length)+" "+str(aAlignment.title.split('|')[4].split(" ")[0])+"\n")
                aFile.write(str(aRecord.sequence)+"\n")
        aFile.close()
        print aFile.name, "saved!"
        theInputFiles.append(aFile.name)
        

def search_distant_homologues(theSequences, theNumIterations):
    for aSequence in theSequences:
        aFilename = "output/"+str(aSequence)+"_fasta.tmp"
        aFile = open(aFilename, "wt")
        aFile.write('>'+str(aSequence)+'\n'+str(theSequences[aSequence])+'\n')
        aFile.close()
        aPsiblastCmd=NcbipsiblastCommandline(db = 'homologues', query = aFilename, evalue =0.001 , out='output/'+str(aSequence)+"_results.xml", outfmt = 5, out_pssm = 'output/'+str(aSequence)+"_pssm.xml", cmd="psiblast", num_iterations=theNumIterations)
        stdout, stderr = aPsiblastCmd()
#       result_handle = open('output/'+str(aSequence)+"_results.xml")
#       brecs = NCBIXML.parse(result_handle)
#       brec = brecs.next()

def align_clustalw(theSequences):
    for aSequence in theSequences:
        cline = ClustalwCommandline("clustalw", infile="output/"+str(aSequence)+"_homologues.fasta")
        stdout, stderr = cline()

################################################################################        
##### PROSITE                                                             ######
################################################################################
def get_prosite_patterns(theAcessionNumber):
    aGoList = []
    handle = ExPASy.get_sprot_raw(theAcessionNumber)
    try:
        record = SwissProt.read(handle)
    except ValueError:
        return "None"
#       print "WARNING: Accession %s not found" % theAcessionNumber
    aList = []
    for aRef in record.cross_references:
        if aRef[0]=='PROSITE':
            aList.append(aRef[1])
    aReturn = str(theAcessionNumber)
    for i in aList:
        aReturn += ' '+str(i)
    return aReturn
    
def search_interproscan_families(theProteins):
    for aSubFamily in theProteins:
        for aProtein in theProteins[aSubFamily]:
            print str(aSubFamily)+' '+get_prosite_patterns(aProtein['uniprot'])

def search_interproscan_subfamilies(theProteins):
    for aProtein in theProteins:
        aSrc = open('output/'+aProtein['uniprot']+'_homologues.fasta', "rt")
        aDst = open('output/'+aProtein['uniprot']+'_homologues_patterns.txt', "wt")
        for aSequence in SeqIO.parse(aSrc,'fasta'):
            aId=aSequence.id.split(' ')[0]
            aDst.write(str(get_prosite_patterns(aId))+'\n')
        aSrc.close()
        aDst.close()

################################################################################        
##### PDB                                                                 ######
################################################################################    
def getPdbId(theProtein):
    aRecord = get_record(theProtein)
    for aRef in aRecord.cross_references:
        if 'PDB' in aRef[0] and 'PDBsum' not in aRef[0]:
            return aRef[1]
    return 'NULL'

def getPdb(theProtein):
    print "\n    BioPython Webservices"
    
    pdb_id = getPdbId(theProtein)
    if 'NULL' not in pdb_id:
        url = 'http://www.rcsb.org/pdb/files/'+pdb_id+'.pdb'
        aFile = open('output/'+pdb_id+'.pdb', "wt")
        print 'Saving the file '+pdb_id+'.pdb'
        for aLine in urllib.urlopen(url).readlines():
            aFile.write(aLine)
        aFile.close()
        print aFile.name, "saved!"  

    
################################################################################
##### Gene Ontology                                                       ######
################################################################################
def get_go_list(theAcessionNumber):
    aQG = QuickGO(verbose=False)
    aGoList = []
    aQueryResult = aQG.Annotation(protein=theAcessionNumber,col="goID,goName", frmt='tsv')
    for aLine in aQueryResult.splitlines():
        teste= aLine.split('\t')
        if 'ID' not in teste[0]:
                #print teste[0], teste[1] # GO:0006629 lipid metabolic process
                aGoList.append((teste[0],teste[1]))
    return aGoList

def print_table_row_subfamily(theProteins,theUniIDs):
    aStr=""
    soma = 0
    aProteins = theUniIDs # ["P32246", "P41597", "P51677", "P51679", "P51681", "P51685", "O00421", "P46094", "P49238", "O60478"]    
    for aProtein in aProteins:
        aStr = aStr+"\t"
        
        if aProtein in theProteins:
            aStr = aStr+"1"
            soma = soma + 1
        else:
            aStr = aStr+" "
    aStr = aStr +"\t"+soma
    return aStr

# Create a CSV file containing a table that correlates GO identifiers with all proteins from a list
def save_table_go_subfamily(theProteins, theFilename, theUniIDs):
    aGoList = {}
    theFilename2 = str(theFilename[:-4])+"_counts.txt"
    aFile = open(theFilename, "wt")
    for aProtein in theProteins:
        #print aProtein
        aGoList[aProtein['uniprot']] = get_go_list(aProtein['uniprot'])
    aDict={}
    for aProtein in aGoList:
        for aGoTuple in aGoList[aProtein]:
            if aGoTuple[0] in aDict:
                aList = aDict[aGoTuple[0]][1]
                aList.append(aProtein)
                aDict[aGoTuple[0]] = (aGoTuple[1], aList)           
            else:
                aDict[aGoTuple[0]] = ((aGoTuple[1], [aProtein]))

    ids = "\t".join(theUniIDs)
    print ids
    aFile.write("GO identifier\tFunction\t"+str(ids)+"\tSUM\n")
    for i in aDict:
        aFile.write(i+"\t"+str(aDict[i][0])+print_table_row_subfamily(aDict[i][1],theUniIDs)+'\n')
    countGOhits(theFilename,theFilename2 ) # sorte by # hits
    aFile.close()
    print aFile.name, " saved!"
    theInputFiles.append(aFile.name)

from collections import *


def countGOhits(theFilename,theFilename2):

        bfile= open(theFilename2, "w")
        aCounter = 0
        new_lines = []
        f = open(theFilename, "r")
        for line in f.readlines()[:1]:
                #print line
                #line = line.strip()
                #items = line.split("\t")
                #for i in items:
                #print i
                #fields = len(items)
                #print line
                #items = items.extend("Counts")
                new_line = "Counts" + "\t" + line
                #new_line.append("\tCounts")
                bfile.write(new_line+"\n")
                #li.insert(2, "new") 

        with open(theFilename, "r") as f:
                for line in f.readlines()[1:]:
                        line = line.strip()
                        items = line.split("\t")
                        
                        for i in items:
                                #print i
                                if i == "1" or "X":
                                        
                                        aCounter += 1

                        #print len(items)#line
                        #fields2 = len(items)
                        #print items[0]    
                        #adj = (fields - fields2)*" \t"
                        #adj = "\t"+aCounter
                        #new_line = "\t".join(line, aCounter)

                        # li.insert(0, 'wow!')
                        #data = file.readlines()
                        #new_line = ( aCounter + line.rstrip('\n') + aCounter + "\n" ) 

                        new_line = str(aCounter)+ "\t" + line 
                        new_lines.append(new_line)
                        
                        bfile.write(new_line+"\n")
                        
                        #file_lines = ['\t'.join([x.strip(), aCounter, '\n']) for x in f.readlines()]
                        aCounter = 0

        #with open(theFilename, 'w') as f:
                #f.writelines(file_lines)
        #new_lines.sort(key=lambda x: x.count, reverse=True)

        num_dict = defaultdict(list)
        for i in new_lines:
                #print i
                items = i.split("\t")
                num = items[0]
                num_dict[num] = items[1:]

        #for i in num_dict.keys():
                #print i, num_dict[i]
                #bfile.write(i+str(num_dict[i])+"\n")

        f.close()
        bfile.close()


def cleanGOtable(theFilename2):#,theFilename3):

        #from CSV import *
        
        count = 0
        count_not = 0
        #records = csv.reader(open(theFilename2), delimiter="\t")
        #for i in records[2]:
                #print i
                
        columns = []
        column = []
        with open(theFilename2) as b:
                for line in b.readlines():
                        if len(line) > 1:
                                fields = line.strip().split("\t")
                                print fields
                                
                                
                        else:
                                line = line.strip()
                                #print line
                                #items = line.split("\t")
                                #columns = line.split("\t")
                                #print columns[2]

                                #col = 2 # third column
                                #filename = '4columns.txt'
                                for col in range(len(fields)):
                                        try:
                                                item = [line.split('\t')[col]] #for line in open(filename,'r')]
                                                column.append(item)
                                                print column
                                                if "X" in columns[1:]:
                                                        print columns
                                                        print i, columns[i], " got it"# GO ID
                                                        count += 1
                                                        next
                                                else:
                                                        print "looking for X"
                                        except:
                                                print "Error: "
                                                
                        
                        #for i in xrange(len(columns)): # reads the column by index
                         #       if "X" in columns[i]:
                          #3              print i, columns[i], " got it"# GO ID
                            #            count += 1
                                        #next
                                #elif "X" not in columns[i]:
                                        #print i, "no GO
                                        #else:
                                                  #print "try later"
                                          #        print "no go", column
                                         #         count_not += 1
                                        
                                
                                        
                                        
                                        
                        #print columns[1]
        print count, count_not

        c += Counter()                  # remove zero and negative counts
                        
                        
def only_gos(theFilename2):
        
        count = 0
        count_not = 0
                
        columns = []
        column = []
        
        col = 2 # third column
        #filename = '4columns.txt'
        third_column = [line[:-1].split('\t')[col] for line in open(theFilename2,'r')]
        print third_column
        


#col_totals = [ sum(x) for x in zip(*my_list) ]
#However, if you need to keep the list,
#[x + y for x, y in zip(*my_list)] #is the fastest.

#for (i, each_line) in enumerate(open('input_file.txt','rb')):
#try:
#column_3 = each_line.split('\t')[2].strip()
#except IndexError:
#print 'Not enough columns on line %i of file.' % (i+1)
#continue

#do_something_with_column_3()
        
def get_column(theFilename2):

        import re
        infile = open ('file', 'r')
        outfile = open('output', 'w')
        column = 31
        for line in infile:
                if not re.match('#', line):
                        line = line.strip()
                        sline = line.split()
                        outfile.write(sline[column] + '\n')
        infile.close()
        outfile.close()
        print outfile.name, " saved!"
        


#countGOhits("1_table_go.txt","1_table_go_counts.txt")

#cleanGOtable("1_table_go_counts.txt")




def countProteinsContainGO(theDict, theGO):
    aCounter = 0
    for aProtein in theDict:
        aFound = False
        for aGoTuple in theDict[aProtein]:
            if theGO in aGoTuple[0] and not aFound:
                aCounter += 1
                aFound = True
    return aCounter

def print_table_row_family(theFamilies):
    aStr=""
    aFamilies = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19'] 
    for aFamily in aFamilies:
        aStr = aStr+"\t"+str(theFamilies[aFamily])
    return aStr

def save_table_go_family(theProteins, theFilename):
#   aDictProtGO={   '1':{ 'CCR1':[GO:03242,GO:38213,GO:329842], 'CCR2':[GO:03242,GO:38213,GO:329842]},
#                   '2':{ 'CCRL1':[GO:45342,GO:435343,GO:3297876], 'CCRL2':[GO:234,GO:33241,GO:343532]}    
#               }
    aDictProtGO={}
    aFile = open(theFilename, "wt")
    print "Getting the GO list of each protein"
    for aSubfamily in theProteins:
        aSubfamilyProteins = theProteins[aSubfamily]
        aGoList = {}
        for aProtein in aSubfamilyProteins:
            aGoList[aProtein['uniprot']] = get_go_list(aProtein['uniprot'])
        aDictProtGO[aSubfamily] = aGoList

#   aDictProtGO={   'GO:32432':{ '1':3, '2': 4 },
#                   'GO:45674':{ '1':5 , '2': 0 }
#               }
    aDictGO = {}
    print "Converting the protein-based list to GO-based dictionary"
    for aSubfamily in aDictProtGO:
        aSubfamilyProteins = aDictProtGO[aSubfamily]
        for aProtein in aSubfamilyProteins:
            for aGoTuple in aSubfamilyProteins[aProtein]:
                if aGoTuple[0] not in aDictGO:
                    aDict = {}
                    for aSubfamily2 in theProteins:
                        aDict[aSubfamily2] = countProteinsContainGO(aDictProtGO[aSubfamily2],aGoTuple[0])
                    aDictGO[aGoTuple[0]] = (aGoTuple[1],aDict)
    print "Writing the CSV file"
    aFile.write("GO identifier\tFunction\tA1\tA2\tA3\tA4\tA5\tA6\tA7\tA8\tA9\tA10\tA11\tA12\tA13\tA14\tA15\tA16\tA17\tA18\tA19"+'\n')
    for i in aDictGO:
        aFile.write(i+"\t"+str(aDictGO[i][0])+print_table_row_family(aDictGO[i][1])+'\n')
    aFile.close()
    print aFile.name, " saved!"
    
    
    
def get_matrix(theFilename, theUniIDs):
    aProteins = theUniIDs # ["CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR8", "CCRL2", "XCR1", "CX3CR1", "GPR137B"]
    aFileIn = open(theFilename, "rt")
    aFileOut = open(str(theFilename[:-4])+"_matrix.txt", "wt")
    aDict = {}
    for aProtein in aProteins:
        aDict[aProtein] = {}
        for aOther in aProteins:
            aDict[aProtein][aOther]=0
    
    for aLine in aFileIn.readlines():
        aRow = aLine.strip().split("\t")
        if "GO:" in aRow[0]:
            aCol=0
            aHitsList = []
            for aElement in aRow[2:12]:
                if 'X' in aElement:
                    aHitsList.append(aProteins[aCol])
                aCol += 1
            for aProt in aHitsList:
                for aOther in aHitsList:
                    if aProt in aOther:
                        aDict[aProt][aOther] += 1
                    else:
                        aDict[aProt][aOther] += 1
    
    ids = "\t".join(theUniIDs)
    aFileOut.write('Protein\t'+str(ids)+'\n')
    for aProtein in aProteins:
        aStr= aProtein+'\t'
        for aOther in aProteins:
            aStr += str(aDict[aProtein][aOther])+'\t'
        aFileOut.write(aStr+'\n')
    aFileIn.close()
    aFileOut.close()
    print aFileOut.name, " saved!"
    theInputFiles.append(aFileOut.name)


def get_hierarchical_candidates(theProteins):
    aGoList = {}
    aCandidates = {}
    aFileOut = open("output/hirarchical_candidates.csv", "wt")
    for aProtein in theProteins:
        aGoList[aProtein['uniprot']] = get_go_list(aProtein['uniprot'])
        aQG = QuickGO(verbose=False)
        aCandidateList = []
        for aGoTerm in aGoList[aProtein['uniprot']]:
            aQueryResult = aQG.Term(aGoTerm[0], frmt='obo')
            for aLine in aQueryResult.splitlines():
                teste= aLine.split(' ')
                if 'is_a:' in teste[0]:
                    aCandidateList.append(teste[1])
        aCandidates[aProtein['uniprot']] = aCandidateList
    for aProtein in theProteins:
        aFileOut.write('##### '+str(aProtein['uniprot'])+'\n')
        for aGoTerm in aCandidates[aProtein['uniprot']]:
            if aGoTerm not in aGoList[aProtein['uniprot']]:
                aFileOut.write(aGoTerm+', '+'\n')

    print aFileOut.name, " saved!"
    theInputFiles.append(aFileOut.name)
    
        
            
            
################################################################################
##### PubMed                                                              ######
################################################################################
def save_references(theAccessions):
    for aAccession in theAccessions:
        aFile = open("output/"+str(aAccession['name'])+"_references.txt", "wt")
        aRecord = get_record(aAccession['uniprot'])
        for aReference in aRecord.references:
            for aRegister in aReference.references:
                if "PubMed" in aRegister[0]:
                    aHandle = Entrez.efetch(db="pubmed", id=aRegister[1], rettype="medline",retmode="text")
                    aArticles = Medline.parse(aHandle)
                    for aArticle in aArticles:
                        if aArticle is not None:
                            if aArticle.has_key('AB'):
                                aFile.write("ID: "+str(aArticle["PMID"])+'\n')
                                aFile.write("TI: "+str(aArticle["TI"])+'\n')
                                aFile.write("AU: "+str(aArticle["AU"])+'\n')
                                aFile.write("AB: "+str(aArticle["AB"])+'\n')
                                aFile.write('\n')
        aFile.close()
        print aFile.name, " saved!"
        theInputFiles.append(aFile.name)
        print "\nBiopython - PDB - ACCESS PubMed - DONE!\n"
        
###########################################################################
#############################################################################3
# 25 mar 2015

#######################

def get_bioblast_record(Unique_Accessions):


    time_begin = datetime.datetime.fromtimestamp(time.time())

    global theInputFiles

    print "\nSearch Accession in Database : SwissProt with Biopython"
        
    
    aFile = open("output/"+str(Unique_Accessions[0])+"_bp_swissprot.txt", "wt")
    header = "Entry\tEntry name\tInterPro\tGene ontology IDs\tec\tKegg Ortholog\tOrganism\tOrganism ID\n"
    aFile.write(header)
    print "\nTotal Accessions : ", Unique_Accessions
    print "Retrieving data in format :\n %s" % header
    iprs = []
    gos = []

    ipprs = []
    goos = []
    
    
    for aProtein in Unique_Accessions:# len = 5355
        handle = ExPASy.get_sprot_raw(aProtein)
        seq_record = SeqIO.read( handle,"swiss")
        handle.close()
        #print dir(seq_record)
        #print "\n\n", seq_record.id 
        ac = seq_record.id 
        #print seq_record.name
        #print seq_record.description
        desc =  seq_record.description
        ec = ""
        ko = ""
        for i in desc.split("; "):
                #print i
                if i.startswith("RecName: "):
                        names = i.split(";")[-1]
                        #print names
                        name = names.split("=")[-1]
                        #print name
                elif i.startswith("SubName: "):
                        names = i.split("=")[-1]
                        name = names.split("{")[0]
                        #print name
                elif i.startswith("EC="):
                        ec = i
            
        for aref in seq_record.dbxrefs:
            #print aref
            if aref[0:3] == "GO:":
                #print aProtein, aref[3:]
                go = aref[3:]
                gos.append(go)
                # print "authors:", aref.authors deprecated
                # print "title:", aref.title
            elif aref[0:9] == "InterPro:":
                ipr = aref[9:]
                iprs.append(ipr)
            elif aref[0:3] == "KO:":
                ko = aref[3:]
                      
        #print seq_record.seq
        #print "Length %i" % len(seq_record)
        """ comment, date_last_sequence_update, taxonomy, date_last_annotation_update, keywords, references,
        accessions, ncbi_taxid, date, organism, gene_name
            """
        #for i in seq_record.annotations["references"]: print i
        #for i in seq_record.annotations["keywords"]: print i

        
            
        org_name =  seq_record.annotations["organism"]
        for i in seq_record.annotations["ncbi_taxid"]:
            tax_id = i

        #print seq_record.annotations["taxonomy"] # lineage
        if ec == "":
                ec = "null"
        if ko == "":
                ko = "null"
        aFile.write(ac+"\t"+name+"\t"+"; ".join(iprs)+"\t"+"; ".join(gos)+"\t"+ec+"\t"+ko+"\t"+org_name+"\t"+tax_id+"\n")

        ipprs = ipprs + iprs
        goos = goos + gos
        del iprs[:]
        del gos[:]

    #c  = Counter(accs)
    #print "\nTOTAL ACS found : ", len(accs) # 6597
    #print "UNIQUE ACS : ", len(c), "\n", c.most_common(10) # 2088
    #for i in dict(c.most_common(10)):
     #   print i, c[i]

     #zip(m.values(), m.keys())
    #[(1, 'a'), (3, 'c'), (2, 'b'), (4, 'd')]
    #>>> mi = dict(zip(m.values(), m.keys()))

    d  = Counter(ipprs)
    print "\nTOTAL IPRs found : ", len(ipprs) # 6597
    print "UNIQUE IPRs : ", len(d), "\n", d.most_common(10) # 3606
    
    #for i in d.most_common(25):
        #print i

    e  = Counter(goos)
    print "\nTOTAL IPRs NAMES found : ", len(goos) # 6597
    print "UNIQUE IPR NAMES : ", len(e), "\n", e.most_common(10) # 3606
    #for i in e.most_common(25):
        #print i
        

    aFile.close()
    #global theInputFiles
    theInputFiles.append(aFile.name)
    print "\n", aFile.name, " saved !"
    print "\nBioBlast - ExPASy - ACCESS ACC - DONE!\n\n"


    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    #winsound.Beep(Freq,Dur)
    print "(complete!)"
    time_end = datetime.datetime.fromtimestamp(time.time())
    print("Time elapsed: ", str(time_end - time_begin))
        

#ids = ['E9AHJ2', 'A4HRI5', 'A4I021', 'A4I9K9', 'A4IE50', 'A4HRQ0', 'A4HY35', 'A4HZ45', 'A4I1E1']
#get_bioblast_record(ids)


##############################################################################################
# BIOPYTHON
##########################################################################################

def get_record1(theAcessionNumber):
	handle = ExPASy.get_sprot_raw(theAcessionNumber)
	try:
		return SwissProt.read(handle)
	except ValueException:
		print "WARNING: Accession %s not found" % theAcessionNumber



        
def getPdbinfo(theProteins):

        from collections import defaultdict

        # name = raw_input("Enter file name : ")
        name = random.choice(theProteins)
        f1 = open("output/"+ str(name)+"_pdb_info.txt", "w")
        f1.write("UniProt ID\tReferences\tID\tID Name\n")
        the_Refs = []
        components = defaultdict(list)
        for i in theProteins:
                aRecord = get_record1(i)
                #print "hello:: Got PDBid"
                for aRef in aRecord.cross_references:
                        #print i, aRef
                        crefs = list(aRef)
                        #print crefs
                        #aRefs.append(aRef)
                        if "InterPro" in aRef[0]:
                                f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        elif "GO" in aRef[0]:
                                f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")  #str(aRef[:])
                                cat = aRef[2]
                                if cat[0] == "C":
                                    components[i].append(aRef[2])
                                    print i, aRef[2]
                            
                        elif 'Proteomes' in aRef[0]:
                                #the_Refs.append()
                                f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        elif "EC" in aRef[0]:
                                f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                                                        
                        #elif "KO" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        elif "BRENDA" in aRef[0]:
                                f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        
    
        
                        
        #the_Refs = " |".join(aRefs)
        f1.close()
        
        #for i , j in components.iteritems():
        #print i, j
       
        print f1 , "is saved!"
        theInputFiles.append(f1.name)
        print "\nBioservices - PDB - ACCESS PDB INFO- DONE!\n"
        
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print "(complete!)"
        #return components
        

#uni_ids= ["A4HT44","A4HUG0","A4IDF2","A4I5K8","A4I066","A4I067"]
#getPdbinfo(uni_ids)


##################################################################################################
##### SWISS PROT - GO ANNOTATIONS
#####                               CELL  COMPOENT - BIOLOGICAL PROCESS - MOLECULAR FUNCTION 
##################################################################################################


def getGOcomponent(filename):

        from collections import defaultdict
        global theInputFiles

        # name = raw_input("Enter file name : ")
        letters= ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
        numbers = ["1","2","3","4","5","6","7","8","9","0"]

        print "Choice =" , filename
        theProteins = []
        f0 = open(filename,"r")

        

        count = 0
        header = f0.readline()
        print "Header\n",header
        #print header
        for i in f0.readlines():
            i = i.strip()
            items = i.split("\t")
            line = items[0]
            if line[0] in letters and line[1] in numbers and line[-1] in numbers and line[-2]: # A4HUD2
                #for i in items:
                #ac = items[0]
                if line not in theProteins:
                    theProteins.append(line)
                    count += 1
            else:
                pass
          

        ###    EXPECTED   TIME
        alph = 0.003 # 4000 proteins / 720 secs
        expected_time = count * float(alph)
        expected_time = "{0:.1f}".format(expected_time)
        print "TOTAL PROTEINS :", len(theProteins)
        print "Expected Time : ", expected_time , "minutes "
        ###
        
        
        name2 = str(filename[:-4])
        #name = random.choice(theProteins)
        f1 = open(name2+"_go_terms_cc.txt", "w")
        f1.write("UniProt ID\tReferences\tID\tID Name\tSource\n")

        
        f2 = open(name2+"_go_component.txt", "w")
        f2.write("UniProt ID\tReferences\tID\tID Name\tSource\n")
        


        the_Refs = []
        components = defaultdict(list)
        gocomp_dict = defaultdict(list)
        all_components = []
        all_go_components = []
        for i in theProteins:
                #print i
                aRecord = get_record1(i) # SWISSPROT - BIOPYTHON
                #print "hello:: Got PDBid"
                for aRef in aRecord.cross_references:
                        #print i, aRef
                        crefs = list(aRef)
                        #print crefs
                        #aRefs.append(aRef)
                        #if "InterPro" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        if "GO" in aRef[0]:
                                f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")  #str(aRef[:])
                                cat = aRef[2]
                                goid = aRef[1]
                                if cat[0] == "C":
                                    components[i].append(aRef[2])
                                    all_components.append(aRef[2])
                                    all_go_components.append(aRef[1])
                                    gocomp_dict[goid].append(aRef[2])
                                    #print i, aRef[2]
                                    f2.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                                
                                
                            
                        #elif 'Proteomes' in aRef[0]:
                        #the_Refs.append()
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        #elif "EC" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                                                        
                        #elif "KO" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        #elif "BRENDA" in aRef[0]:
                                #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        
    
        
                        
        #the_Refs = " |".join(aRefs)
        f1.close()
        f2.close()
        
        #for i , j in components.iteritems():
        #print i, j
        print "all proteins : ", len(theProteins), "\nkeys : ", len(components.keys())

        # PIE PLOT #
        
        import matplotlib.pyplot as plt
        from collections import Counter
        from collections import OrderedDict
        # http://matplotlib.org/examples/pie_and_polar_charts/pie_demo_features.html
        
        cnt = Counter(all_components)
        cnt2 = Counter(all_go_components)
        print cnt.most_common(25)

        # The slices will be ordered and plotted counter-clockwise.
        labels = []
        sizes = []
        f3 = open(str(filename[:-4])+"_sum_components.txt","w")

        mydict = OrderedDict(sorted(cnt.most_common(), key=lambda t: t[1]),reverse=True)
        mydict2 = OrderedDict(sorted(cnt2.most_common(), key=lambda t: t[1]),reserse=True)
        
        
        for i,j in mydict2.iteritems():# cnt.most_common():
            f3.write(str(i)+"\t"+str(gocomp_dict[i])+"\t"+str(j)+"\n")
            #print i, j
            
        for i,j in cnt2.most_common(25):
            #labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
            #sizes = [15, 30, 45, 10]
            labels.append(i[:50])
            sizes.append(j)

        label_hits = []
        for i,j in cnt.most_common(25):# cnt.most_common():
            nam = str(j) +" "+ str(i[:50])
            label_hits.append(nam)
            
        tot_acs = len(theProteins)
        attr_acs = len(components.keys())
        left_acs = tot_acs - attr_acs
        print "labels, attr_acs, tot_acs"
        print len(labels), attr_acs, tot_acs
        #labels.append("Unknown")
        #sizes.append(left_acs)
        print "Unknown : ", left_acs, "\t", float(left_acs)/tot_acs*100, "%"
        colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
        explode = (0, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

        fig = plt.figure()
        graph_title = str(filename[:-4]) + "_"+ str(len(components.keys()))
        fig.canvas.set_window_title(graph_title)

        #lg = legend()
        #lg.draw_frame(False)
        #ax.legend([extra, bar_0_10, bar_10_100], ("My explanatory text", "0-10", "10-100"))

        legg = str("Acs : "+len(theProteins)+ " keys : "+ len(components.keys()))
        plt.suptitle("Top 25 - GO Cell Component"+legg)

        plt.pie(sizes,  labels=labels, colors=colors,
                        autopct='%1.1f%%', shadow=True, startangle=90) #explode=explode,
        # Set aspect ratio to be equal so that pie is drawn as a circle.
        plt.axis('equal')

        plt.legend(label_hits, loc='best', shadow=True)

        print f1.name , " saved!"
        print f2.name, " saved!"
        print f3.name, " saved!"
        theInputFiles.append(f1.name)
        theInputFiles.append(f2.name)
        theInputFiles.append(f3.name)
        print "\nBioservices - SwissProt - ACCESS GO Cellular Component - DONE!\n"
        
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print "(complete!)"
        #return components

        plt.show()
       
        

#getGOcomponent("output/2peps_all_.txt")
#getGOcomponent("output/610_all.txt")
#getGOcomponent("20150831_estes/2peps_3999_filter_207_A4HTV2.txt")


#####       GET GO BIOLOGICAL PROCESS       ###########################################3

def getGOprocess(filename):

        from collections import defaultdict
        global theInputFiles

        # name = raw_input("Enter file name : ")
        letters= ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
        numbers = ["1","2","3","4","5","6","7","8","9","0"]
        

        print "Choice =" , filename
        theProteins = []
        count = 0
        f0 = open(filename,"r")

        header = f0.readline()
        for i in f0.readlines():
            i = i.strip()
            items = i.split("\t")
            line = items[0]
            if line[0] in letters and line[1] in numbers and line[-1] in numbers and line[-2]: # A4HUD2
                #for i in items:
                #ac = items[0]
                if line not in theProteins:
                    theProteins.append(line)
                    count += 1
            else:
                pass
            

        ###    EXPECTED   TIME
        alph = 0.003 # 4000 proteins / 720 secs
        expected_time = count * float(alph)
        expected_time = "{0:.1f}".format(expected_time)
        print "TOTAL PROTEINS :", len(theProteins)
        print "Expected Time : ", expected_time , "minutes "
        ###

        
       
        name2 = str(filename[:-4])
        #name = random.choice(theProteins)
        f1 = open(name2+"_go_terms_bp.txt", "w")
        f1.write("UniProt ID\tReferences\tID\tID Name\tSource\n")

        
        f2 = open(name2+"_go_process.txt", "w")
        f2.write("UniProt ID\tReferences\tID\tID Name\tSource\n")
        
        the_Refs = []
        components = defaultdict(list)
        gocomp_dict = defaultdict(list)
        all_components = []
        all_go_components= []
        for i in theProteins:
                #print i
                aRecord = get_record1(i)
                #print "hello:: Got PDBid"
                for aRef in aRecord.cross_references:
                        #print i, aRef
                        crefs = list(aRef)
                        #print crefs
                        #aRefs.append(aRef)
                        #if "InterPro" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        if "GO" in aRef[0]:
                                f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")  #str(aRef[:])
                                cat = aRef[2]
                                goid = aRef[1]
                                if cat[0] == "P":
                                    components[i].append(aRef[2])
                                    all_components.append(aRef[2])
                                    all_go_components.append(aRef[1])
                                    gocomp_dict[goid].append(aRef[2])
                                    #print i, aRef[2]
                                    f2.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                                    
                                
                                
                            
                        #elif 'Proteomes' in aRef[0]:
                        #the_Refs.append()
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        #elif "EC" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                                                        
                        #elif "KO" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        #elif "BRENDA" in aRef[0]:
                                #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        
    
        
                        
        #the_Refs = " |".join(aRefs)
        f1.close()
        f2.close()
        
        #for i , j in components.iteritems():
        #print i, j
        print "all proteins : ", len(theProteins), "\nkeys : ", len(components.keys())

        # PIE PLOT #
        
        import matplotlib.pyplot as plt
        from collections import Counter
        from collections import OrderedDict
        # http://matplotlib.org/examples/pie_and_polar_charts/pie_demo_features.html
        
        cnt = Counter(all_components)
        print cnt.most_common(25)
        cnt2 = Counter(all_go_components)

        # The slices will be ordered and plotted counter-clockwise.
        labels = []
        sizes = []
        f3 = open(str(filename[:-4])+"_sum_process.txt","w")

        mydict = OrderedDict(sorted(cnt.most_common(), key=lambda t: t[1]))
        mydict2 = OrderedDict(sorted(cnt2.most_common(), key=lambda t: t[1]))
        
        for i,j in mydict2.iteritems():# cnt.most_common():
            f3.write(str(i)+"\t"+str(gocomp_dict[i])+"\t"+str(j)+"\n")
            #f3.write(str(i)+"\t"+str(j)+"\n")
            #print i, j
            
        for i,j in cnt2.most_common(25):
            #labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
            #sizes = [15, 30, 45, 10]
            labels.append(i[:50])
            sizes.append(j)
            
        label_hits = []
        for i,j in cnt.most_common(25):# cnt.most_common():
            nam = str(j) +" "+ str(i[:50])
            label_hits.append(nam)
            
        tot_acs = len(theProteins)
        attr_acs = len(components.keys())
        left_acs = tot_acs - attr_acs
        print "labels, attr_acs, tot_acs"
        print len(labels), attr_acs, tot_acs
        #labels.append("Unknown")
        #sizes.append(left_acs)
        print "Unknown : ", left_acs, "\t", float(left_acs)/tot_acs*100, "%"
        colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
        explode = (0, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

        fig = plt.figure()
        graph_title = str(filename[:-4]) + "_"+ str(len(components.keys()))
        fig.canvas.set_window_title(graph_title)

        #lg = legend()
        #lg.draw_frame(False)
        #ax.legend([extra, bar_0_10, bar_10_100], ("My explanatory text", "0-10", "10-100"))

        legg = str("Acs : "+len(theProteins)+ " keys : "+ str(len(components.keys())))
        plt.suptitle("Top 25 - GO Biological Process"+legg)

        plt.pie(sizes,  labels=labels, colors=colors,
                        autopct='%1.1f%%', shadow=True, startangle=90) #explode=explode,
        # Set aspect ratio to be equal so that pie is drawn as a circle.
        plt.axis('equal')

        plt.legend(label_hits, loc='best', shadow=True)

        print f1.name , " saved!"
        print f2.name, " saved !"
        print f3.name, " saved !"
        theInputFiles.append(f1.name)
        theInputFiles.append(f2.name)
        theInputFiles.append(f3.name)
        print "\nBioservices - SwissProt - ACCESS GO Biological Process - DONE!\n"
        
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print "(complete!)"
        #return components

        plt.show()
       
        

#getGOprocess("20150831_estes/2peps_3999_filter_207_A4HTV2.txt")


#####       GET GO MOLECULAR FUNCTION       ###########################################3

def getGOfunction(filename):

        from collections import defaultdict
        global theInputFiles

        # name = raw_input("Enter file name : ")
        letters= ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
        numbers = ["1","2","3","4","5","6","7","8","9","0"]

        print "Choice = " , filename
        theProteins = []
        count = 0
        f0 = open(filename,"r")
        header = f0.readline()
        for i in f0.readlines():
            i = i.strip()
            items = i.split("\t")
            line = items[0]
            if line[0] in letters and line[1] in numbers and line[-1] in numbers and line[-2]: # A4HUD2
                #for i in items:
                #ac = items[0]
                if line not in theProteins:
                    theProteins.append(line)
                    count += 1
            else:
                pass
            #print ac

        ###    EXPECTED   TIME
        alph = 0.003 # 4000 proteins / 720 secs
        expected_time = count * float(alph)
        expected_time = "{0:.1f}".format(expected_time)
        print "TOTAL PROTEINS :", len(theProteins)
        print "Expected Time : ", expected_time , "minutes "
        ###

        name2 = str(filename[:-4])
        #name = random.choice(theProteins)
        f1 = open(name2+"_go_terms_mf.txt", "w")
        f1.write("UniProt ID\tReferences\tID\tID Name\tSource\n")

        
        f2 = open(name2+"_go_function.txt", "w")
        f2.write("UniProt ID\tReferences\tID\tID Name\tSource\n")

        
        the_Refs = []
        components = defaultdict(list)
        gocomp_dict = defaultdict(list)
        all_components = []
        all_go_components = []
        for i in theProteins:
                #print i
                aRecord = get_record1(i)
                #print "hello:: Got PDBid"
                for aRef in aRecord.cross_references:
                        #print i, aRef
                        crefs = list(aRef)
                        #print crefs
                        #aRefs.append(aRef)
                        #if "InterPro" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        if "GO" in aRef[0]:
                                f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")  #str(aRef[:])
                                cat = aRef[2]
                                goid = aRef[1]
                                if cat[0] == "F":
                                    components[i].append(aRef[2])
                                    all_components.append(aRef[2])
                                    all_go_components.append(aRef[1])
                                    gocomp_dict[goid].append(aRef[2])
                                    f2.write(str(i)+"\t"+"\t".join(crefs)+"\n")  #str(aRef[:])
                                    
                                
                                    
                                    #print i, aRef[2]
                                
                                
                            
                        #elif 'Proteomes' in aRef[0]:
                        #the_Refs.append()
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        #elif "EC" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                                                        
                        #elif "KO" in aRef[0]:
                        #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        #elif "BRENDA" in aRef[0]:
                                #f1.write(str(i)+"\t"+"\t".join(crefs)+"\n")
                        
    
        
                        
        #the_Refs = " |".join(aRefs)
        f1.close()
        f2.close()
        
        #for i , j in components.iteritems():
        #print i, j
        print "all proteins : ", len(theProteins), "\nkeys : ", len(components.keys())

        # PIE PLOT #
        
        import matplotlib.pyplot as plt
        from collections import Counter
        from collections import OrderedDict
        # http://matplotlib.org/examples/pie_and_polar_charts/pie_demo_features.html
        
        cnt = Counter(all_components)
        cnt2 = Counter(all_go_components)
        print cnt.most_common(25)

        # The slices will be ordered and plotted counter-clockwise.
        labels = []
        sizes = []
        f3 = open(str(filename[:-4])+"_sum_function.txt","w")

        mydict = OrderedDict(sorted(cnt.most_common(), key=lambda t: t[1]),reverse=True)

        mydict2 = OrderedDict(sorted(cnt2.most_common(), key=lambda t: t[1]),reverse=True)

       
        for i,j in mydict.iteritems():# cnt.most_common():
            f3.write(str(i)+"\t"+str(gocomp_dict[i])+"\t"+str(j)+"\n")
            #f3.write(str(i)+"\t"+str(j)+"\n")
            #print i, j
            
        for i,j in cnt2.most_common(25):
            #labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
            #sizes = [15, 30, 45, 10]
            labels.append(i[:50])
            sizes.append(j)

        label_hits = []
        for i,j in cnt.most_common(25):# cnt.most_common():
            nam = str(j) +" "+ str(i[:50])
            label_hits.append(nam)
            
        tot_acs = len(theProteins)
        attr_acs = len(components.keys())
        left_acs = tot_acs - attr_acs
        print "labels, attr_acs, tot_acs"
        print len(labels), attr_acs, tot_acs
        #labels.append("Unknown")
        #sizes.append(left_acs)
        print "Unknown : ", left_acs, "\t", float(left_acs)/tot_acs*100, "%"
        colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
        explode = (0, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

        fig = plt.figure()
        graph_title = str(filename[:-4]) + "_"+ str(len(components.keys()))
        fig.canvas.set_window_title(graph_title)

        #lg = legend()
        #lg.draw_frame(False)
        #ax.legend([extra, bar_0_10, bar_10_100], ("My explanatory text", "0-10", "10-100"))

        legg = str("Acs : "+len(theProteins)+ " keys : "+ len(components.keys()))
        plt.suptitle("Top 25 - GO Molecular Function"+legg)

        plt.pie(sizes,  labels=labels, colors=colors,
                        autopct='%1.1f%%', shadow=True, startangle=90) #explode=explode,
        # Set aspect ratio to be equal so that pie is drawn as a circle.
        plt.axis('equal')

        plt.legend(label_hits, loc='best', shadow=True)

        print f1.name , " saved!"
        print f2.name, "saved!"
        print f3.name, " saved!"
        theInputFiles.append(f1.name)
        theInputFiles.append(f2.name)
        theInputFiles.append(f3.name)
        print "\nBioservices - SwissProt - ACCESS GO Molecular Function - DONE!\n"
        
        Freq = 1000 # Set Frequency To 2500 Hertz
        Dur = 100 # Set Duration To 1000 ms == 1 second
        winsound.Beep(Freq,Dur)
        print "(complete!)"
        #return components

        plt.show()
       
        

#getGOfunction("20150831_estes/2peps_3999_filter_207_A4HTV2.txt")
