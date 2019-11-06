import datetime, time
time_begin = datetime.datetime.fromtimestamp(time.time())


from bi_biopython_functions_tk import *

def access_blast_with_biopython3(theUniIDs):#theTaxIDs
        # 0) Write UniProt IDs to File
        
        #theFilename = "GPCR_proteins1.txt"
        theFilename = "1_name_ac.txt"

        ##write_ids_to_file(theUniIDs,theFilename)
        from_uni_id_get_bs_uniprot(theUniIDs, theFilename)

        # 1) Load all proteins accession numbers to a dictionary
        print ('[1] Loading protein accessions')
        allProteins = load_proteins_accessions("1_name_ac.txt") # search file
        # 2) Initialize a variable only with the proteins from organism ID (requires 1)
        print ('[2] Selecting protein accessions from Leishmania Infantum')
        LI_Proteins = allProteins[theTaxID]#
        # 3) Retrieve and save the aminoacid sequences of the proteins (requires 2)
        print ('[3] Retrieving and saving proteins sequences')
        save_proteins_sequences(LI_Proteins, "output/LI_proteins_sequence.fasta")
        # 4) Load from the file the amino-acid sequences of the proteins (requires 3)
        print ('[4] Loading proteins sequences')
        LI_Sequences = load_proteins_sequences("output/LI_proteins_sequence.fasta")
        # 5) Retrieve and save the N homologues of each protein sequence (requires 4)
        print ('[run.py] Retrieving and saving homologues proteins')
        save_homologues(LI_Sequences,10)
        # 6) Retrieve and save the sequences of distant homologues of each protein sequence (requires 4)
        print ('[run.py] Retrieving and saving distant homologues proteins')
        # search_distant_homologues(LI_Sequences, 10)
        # 7) Retrieve and save the abstract of references related to each protein (requires 2)
        print ('[run.py] Retrieving and saving abstracts from cross-references')
        #save_references(LI_Proteins)
        # 8) Perform a multiple sequence alignment using the Clustal-W algorithm (requires 4 and 5)
        print ('[run.py] Performing multiple sequence alignments with Clustal-W')
        #align_clustalw(LI_Sequences)
        # 9) Retrieve and print the patterns annotated in InterProScan regarding a family of proteins (requires 1 and 5)
        print ('[run.py] Retrieving patterns from InterProScan (for all families)')
        #search_interproscan_families(allProteins)
        # 10) Retrieve and print the patterns annotated in InterProScan regarding a subfamily of proteins (requires 1 and 5)
        print ('[run.py] Retrieving patterns from InterProScan (for subfamily homologues)')
        ###search_interproscan_subfamilies(LI_Proteins)
        # 11) Retrieve GO terms and create a associative table for a list of proteins (requires 2)
        print ('[run.py] Creating GO table for subfamily A1')
        ###save_table_go_subfamily(LI_Proteins, "output/subfamily_A1_table_go.csv")
        # 12) Create a correlation matrix regarding a GO table created in previous steps (requires 11)
        print ('[run.py] Creating correlation matrix (for subfamily A1)')
        ###get_matrix("output/subfamily_A1_table_go.csv")
        # 13) Retrieve GO terms and create a associative table for a list of subfamilies (requires 1)
        print ('[run.py] Creating GO table for family A')
        ###save_table_go_family(allProteins, "output/family_A_table_go.csv")
        # 14) Retrieve and print new candidate GO terms based in hierarchical analysis (requires 2)
        print ('[run.py] Running hierarchical suggestion algorithms for GO terms')
        #get_hierarchical_candidates(LI_Proteins)
        # 15) Retrieve and save the PDB structure (if exist) of a protein (requires 2)
        print ('[run.py] Retrieving ID and PDF files')
        for i in LI_Proteins:
                getPdbId(str(i['uniprot']))
                getPdb(str(i['uniprot']))

	# 16) Retrieve and save the cross-references from ExPASy SwissProt - BioPython
	###save_cross_references(LI_Proteins,"output/LI_protein_cross_references.txt")
	# 17) Filter cross-references in accession - interpro - go

	#count_ac("output/LI_protein_cross_references.txt")#11143
	#print "done!"
	#print "check output folder. Bye"

#access_blast_with_biopython("5671")
	###

#theUniIDs = ["A4HW10","A4ID37","A4HUD4","A4I2W4"]
#access_blast_with_biopython3(theUniIDs)


def from_uni_ids_get_pdb(uni_ids):
        for i in uni_ids:
                getPdbId(str(i))
                getPdb(str(i))



time_end = datetime.datetime.fromtimestamp(time.time())

print("Time elapsed: ", str(time_end - time_begin))
print ("Thank you!")
