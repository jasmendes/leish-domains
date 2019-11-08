############################################################################
from __future__ import division
import datetime, time
time_begin = datetime.datetime.fromtimestamp(time.time())

import sys, os
from glob import *
#print "FILE PATH:", sys.argv[0]
#print "FILE NAME:", os.path.basename(sys.argv[0])

# /kaggle/input/02-exp.txt



from bioservices import *




from bioservices import QuickGO
from math import pi
import urllib




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

ab = from_tax_id_get_bs_uniprot('5671')
