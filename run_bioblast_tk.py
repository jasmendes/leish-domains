
import sys, os







my_paths1 = "c:\Users\mendes\appdata\local\enthought\canopy\user\scripts" # ;C:\Users\mendes\AppData\Local\Enthought\Canopy\App\appdata\canopy-1.5.2.2785.win-x86_64\DLLs;C:\Users\mendes\AppData\Local\Enthought\Canopy\App\appdata\canopy-1.5.2.2785.win-x86_64\lib;C:\Users\mendes\AppData\Local\Enthought\Canopy\App\appdata\canopy-1.5.2.2785.win-x86_64\lib\plat-win;C:\Users\mendes\AppData\Local\Enthought\Canopy\App\appdata\canopy-1.5.2.2785.win-x86_64\lib\lib-tk"

my_paths2 = "C:\Users\mendes\AppData\Local\Enthought\Canopy\App\appdata\canopy-1.5.2.2785.win-x86_64\lib\site-packages"

sys.path.append(my_paths1)
sys.path.append(my_paths2)
#sys.path.append("C:\Python276")
a, b, c, d = "C:\Python273\Lib\site-packages", "C:\Python276\Lib\site-packages", "C:\Python279\Lib\site-packages","C:\Users\mendes\AppData\Local\Enthought\Canopy\User\lib\site-packages"
sys.path.append(a)
sys.path.append(b)
sys.path.append(c)
sys.path.append(d)

#for i in sys.path:
        #print i

from bi_biopython_functions_tk import *
#import tk_tax_id_v151

from random import randint

def access_bioblast_go_matrix(theUniIDs):
        
        import datetime, time
        time_begin = datetime.datetime.fromtimestamp(time.time())
        
        # 0) Write UniProt IDs to File
        
        #theFilename = "GPCR_proteins1.txt"
        name = random.randint(1,1000)
        theFilename = "output/"+str(name)+"1_name_ac.txt"
        theFilename2 = str(theFilename[:-4])+"_go_matrix.txt"

        ##write_ids_to_file(theUniIDs,theFilename)
        from_uni_id_get_bs_uniprot(theUniIDs, theFilename)
        
        # 1) Load all proteins accession numbers to a dictionary
        print '> Loading protein accessions'
        allProteins = load_proteins_accessions(theFilename) # "GPCR_proteins1.txt")
        # 2) Initialize a variable only with the proteins from the A1 subfamily (requires 1)
        print '> Selecting protein accessions'# from subfamily A1'
        a1Proteins = allProteins['1']
        

        ## 11) Retrieve GO terms and create a associative table for a list of proteins (requires 2)
        print '> Creating GO table 1'
        save_table_go_subfamily(a1Proteins, theFilename2, theUniIDs)
        ## 12) Create a correlation matrix regarding a GO table created in previous steps (requires 11)
        print '> Creating correlation matrix '# (for subfamily A1)'
        get_matrix(theFilename2,theUniIDs)

        ## 13) Retrieve GO terms and create a associative table for a list of subfamilies (requires 1)
        #print '> Creating GO table for family A'
        #save_table_go_family(allProteins, "output/family_A_table_go.csv")
        ## 14) Retrieve and print new candidate GO terms based in hierarchical analysis (requires 2)
        #print '> Running hierarchical suggestion algorithms for GO terms'
        #get_hierarchical_candidates(a1Proteins)

        ## 15) Retrieve and save the PDB structure (if exist) of a protein (requires 2)
        #print '> Retrieving ID and PDB files'
        #for i in a1Proteins:
                #getPdbId(str(i['uniprot']))
                #getPdb(str(i['uniprot']))

        #del theUniIDs[:]
        # RUNS IN 3 Minutes
        time_end = datetime.datetime.fromtimestamp(time.time())
        print("Time elapsed: ", str(time_end - time_begin))

#theUniIDs = ["A4HW10","A4ID37","A4HUD4","A4I2W4"] #,"A4HTB4","A4ICJ6","A4I2Y3","A4IDF2","A4HX89"] # t1 =1min40sec

#access_bioblast_go_matrix(theUniIDs)

####################################################################################

def access_bioblast_abstracts(theUniIDs):
        
        import datetime, time
        time_begin = datetime.datetime.fromtimestamp(time.time())
        
        # 0) Write UniProt IDs to File
        
        #theFilename = "GPCR_proteins1.txt"
        name = random.randint(1,1000)
        theFilename = str(name)+"1_name_ac.txt"

        ##write_ids_to_file(theUniIDs,theFilename)
        from_uni_id_get_bs_uniprot(theUniIDs, theFilename)
        
        # 1) Load all proteins accession numbers to a dictionary
        print '> Loading protein accessions'
        allProteins = load_proteins_accessions(theFilename) # "GPCR_proteins1.txt")
        # 2) Initialize a variable only with the proteins from the A1 subfamily (requires 1)
        print '> Selecting protein accessions'# from subfamily A1'
        a1Proteins = allProteins['1']
        # 3) Retrieve and save the aminoacid sequences of the proteins (requires 2)
        #print '> Retrieving and saving proteins sequences'
        #save_proteins_sequences(a1Proteins, "output/1_proteins_sequence.fasta")
        # 4) Load from the file the aminoacid sequences of the proteins (requires 3)
        #print '> Loading proteins sequences'
        #gSequences = load_proteins_sequences("output/1_proteins_sequence.fasta")
        # 5) Retrieve and save the N homologues of each protein sequence (requires 4)
        #print '> Retrieving and saving homologues proteins'
        #save_homologues(gSequences,10)# takes long time ????

        ## 6) Retrieve and save the sequences of distant homologues of each protein sequence (requires 4)
        #print '> Retrieving and saving distant homologues proteins'
        #search_distant_homologues(gSequences, 10)

        ## 7) Retrieve and save the abstract of references related to each protein (requires 2)
        print '> Retrieving and saving abstracts from cross-references'
        save_references(a1Proteins)

        ## 8) Perform a multiple sequence alignment using the Clustal-W algorithm (requires 4 and 5)
        #print '> Performing multiple sequence alignments with Clustal-W'
        #align_clustalw(gSequences) # ERROR no command

        ## 9) Retrieve and print the patterns annotated in InterProScan regarding a family of proteins (requires 1 and 5)
        #print '> Retrieving patterns from InterProScan (for all families)'
        #search_interproscan_families(allProteins)

        ## 10) Retrieve and print the patterns annotated in InterProScan regarding a subfamily of proteins (requires 1 and 5)
        #print '> Retrieving patterns from InterProScan '#(for A1 subfamily homologues)'
        #search_interproscan_subfamilies(a1Proteins)

        ## 11) Retrieve GO terms and create a associative table for a list of proteins (requires 2)
        #print '> Creating GO table 1'
        #save_table_go_subfamily(a1Proteins, "output/1_table_go.txt", theUniIDs)

        ## 12) Create a correlation matrix regarding a GO table created in previous steps (requires 11)
        #print '> Creating correlation matrix '# (for subfamily A1)'
        #get_matrix("output/1_table_go.txt",theUniIDs)

        ## 13) Retrieve GO terms and create a associative table for a list of subfamilies (requires 1)
        #print '> Creating GO table for family A'
        #save_table_go_family(allProteins, "output/family_A_table_go.csv")
        ## 14) Retrieve and print new candidate GO terms based in hierarchical analysis (requires 2)
        #print '> Running hierarchical suggestion algorithms for GO terms'
        #get_hierarchical_candidates(a1Proteins)

        ## 15) Retrieve and save the PDB structure (if exist) of a protein (requires 2)
        #print '> Retrieving ID and PDB files'
        #for i in a1Proteins:
                #getPdbId(str(i['uniprot']))
                #getPdb(str(i['uniprot']))

        #del theUniIDs[:]
        # RUNS IN 3 Minutes
        time_end = datetime.datetime.fromtimestamp(time.time())
        print("Time elapsed: ", str(time_end - time_begin))

#theUniIDs = ["A4HW10","A4ID37","A4HUD4","A4I2W4"] #,"A4HTB4","A4ICJ6","A4I2Y3","A4IDF2","A4HX89"] # t1 =1min40sec

#access_bioblast_abstracts(theUniIDs)

####################################################################################

def access_bioblast_with_biopython1(theUniIDs):
        
        import datetime, time
        time_begin = datetime.datetime.fromtimestamp(time.time())
        
        # 0) Write UniProt IDs to File
        
        #theFilename = "GPCR_proteins1.txt"
        name = random.randint(1,1000)
        theFilename = str(name)+"1_name_ac.txt"

        ##write_ids_to_file(theUniIDs,theFilename)
        from_uni_id_get_bs_uniprot(theUniIDs, theFilename)
        
        # 1) Load all proteins accession numbers to a dictionary
        print '> Loading protein accessions'
        allProteins = load_proteins_accessions(theFilename) # "GPCR_proteins1.txt")
        # 2) Initialize a variable only with the proteins from the A1 subfamily (requires 1)
        print '> Selecting protein accessions'# from subfamily A1'
        a1Proteins = allProteins['1']
        # 3) Retrieve and save the aminoacid sequences of the proteins (requires 2)
        print '> Retrieving and saving proteins sequences'
        save_proteins_sequences(a1Proteins, "output/1_proteins_sequence.fasta")
        # 4) Load from the file the aminoacid sequences of the proteins (requires 3)
        print '> Loading proteins sequences'
        gSequences = load_proteins_sequences("output/1_proteins_sequence.fasta")
        # 5) Retrieve and save the N homologues of each protein sequence (requires 4)
        print '> Retrieving and saving homologues proteins'
        save_homologues(gSequences,10)# takes long time ????

        ## 6) Retrieve and save the sequences of distant homologues of each protein sequence (requires 4)
        #print '> Retrieving and saving distant homologues proteins'
        #search_distant_homologues(gSequences, 10)

        ## 7) Retrieve and save the abstract of references related to each protein (requires 2)
        #print '> Retrieving and saving abstracts from cross-references'
        #save_references(a1Proteins)

        ## 8) Perform a multiple sequence alignment using the Clustal-W algorithm (requires 4 and 5)
        #print '> Performing multiple sequence alignments with Clustal-W'
        #align_clustalw(gSequences) # ERROR no command
        ## 9) Retrieve and print the patterns annotated in InterProScan regarding a family of proteins (requires 1 and 5)
        #print '> Retrieving patterns from InterProScan (for all families)'
        #search_interproscan_families(allProteins)

        ## 10) Retrieve and print the patterns annotated in InterProScan regarding a subfamily of proteins (requires 1 and 5)
        #print '> Retrieving patterns from InterProScan '#(for A1 subfamily homologues)'
        #search_interproscan_subfamilies(a1Proteins)
        ## 11) Retrieve GO terms and create a associative table for a list of proteins (requires 2)

        #print '> Creating GO table 1'
        #save_table_go_subfamily(a1Proteins, "output/1_table_go.txt", theUniIDs)
        ## 12) Create a correlation matrix regarding a GO table created in previous steps (requires 11)
        #print '> Creating correlation matrix '# (for subfamily A1)'
        #get_matrix("output/1_table_go.txt",theUniIDs)

        ## 13) Retrieve GO terms and create a associative table for a list of subfamilies (requires 1)
        #print '> Creating GO table for family A'
        #save_table_go_family(allProteins, "output/family_A_table_go.csv")
        ## 14) Retrieve and print new candidate GO terms based in hierarchical analysis (requires 2)
        #print '> Running hierarchical suggestion algorithms for GO terms'
        #get_hierarchical_candidates(a1Proteins)

        ## 15) Retrieve and save the PDB structure (if exist) of a protein (requires 2)
        #print '> Retrieving ID and PDB files'
        #for i in a1Proteins:
                #getPdbId(str(i['uniprot']))
                #getPdb(str(i['uniprot']))

        #del theUniIDs[:]
        # RUNS IN 3 Minutes
        time_end = datetime.datetime.fromtimestamp(time.time())
        print("Time elapsed: ", str(time_end - time_begin))

#theUniIDs = ["A4HW10","A4ID37","A4HUD4","A4I2W4"] #,"A4HTB4","A4ICJ6","A4I2Y3","A4IDF2","A4HX89"] # t1 =1min40sec
#access_bioblast_with_biopython1(theUniIDs)


