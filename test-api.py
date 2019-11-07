

### Loading libraries                                                    #####
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


def from_gene_id_get_kegg(gene_id):
    
    letters= ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y","Z"]
    numbers = ["1","2","3","4","5","6","7","8","9","0"]
    
    # ORGANISM : hsa
    # KEGG ID : hsa1234
    # PATHWAY : hsa12345

    k = KEGG(verbose=False)

    print ("KEGG  -  Look for relevante Pathways, using gene . . . ")

    print ("Mapping gene to protein ac . . . ")
    u = UniProt()
    res = u.search(str(gene_id), frmt="tab", columns="genes, id")
    for i in res.split("\n"):
        print (i)
        while len(i.split("\t")) > 1 :
            i = i.split()[1]
            #print i
            if i[0] in letters and i[1] in numbers and i[-1] in numbers:
                print (i)
                uni_id = i

    kegg_ids = []

    print ("Mapping protein ac to kegg ID. . . ")
    # http://pythonhosted.org/bioservices/convertor_tutorial.html
    map1  = u.mapping(fr="ACC+ID", to="KEGG_ID", query=str(uni_id))# P_ENTREZGENEID, "UNIGENE_ID, EMBL,GeneID
    print( map1 )# defaultdict(<type 'list'>, {u'P43403': [u'hsa:7535']})

    for i in map1.values():
        i = str(i)
        kegg_id = i[3:-2]
        _id = i[7:-2]
        print (_id) # LINJ_10_0520
        print (i, " = ", kegg_id) # [u'lif:LINJ_10_0520'] = lif:LINJ_10_0520
        k_organism = i[3:6] # 'lif'
        print (k_organism)
        kegg_ids.append(kegg_id)
    print (k.get(str(kegg_id)))
    get_kegg = k.get(str(kegg_id))


    print ("Find  protein . . . ")
    #pname = raw_input("Please enter protein name : ")
    print (k.find(str(k_organism), uni_id))
    get_pname = k.find(str(k_organism), uni_id)

    print ("Find  pathways for this organism . . . ")
    print (k.list("pathway", organism= k_organism))
    res_paths = k.list("pathway", organism= k_organism)
    pathways = [x.split()[0] for x in res_paths.strip().split("\n")]
    print (len(pathways) ,"pathways found in Kegg organism ", k_organism) # as of Nov 2014 -> 94)
    
   
    print (uni_id, "\n",k.find(str(k_organism), str(kegg_id))) # lif:LINJ_14_0700	putative fatty acid elongase (EC:2.3.1.119)

    get_pname = k.find(str(k_organism), str(kegg_id))

    print ("Find pathway from gene . . . ")
    get_path = k.get_pathway_by_gene(str(_id), str(k_organism))
    
    if get_path.__repr__() == "None":
        print (get_path)
        kegg_path = "None"
    else:
        for i in get_path.keys():
            print (i)
            i = str(i)
            if i.isdigit() == True and len(i) == 5:
                print (i)
                kegg_path = i
                kegg_path = k_organism + kegg_path
                print ("\n",uni_id, " = ",kegg_id, "\nPathway found: ", kegg_path)

    print ("Showing Pathway Marker in internet page . . .")
    #print k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
    
    bf = open(str(gene_id)+"_path.txt","w")
    bf.write(res_paths+"\n")
    bf.write(get_pname+"\n")
    bf.write(get_kegg)

    if get_path.__repr__() != "None":
        for i in get_path:
            bf.write(str(i))
    
    print ("Save image of the pathway . . . ")
    # NF-kappa B signaling pathway - Homo sapiens (human)
    # res0 =  k.get("hsa05130/image")
    if kegg_path != "None":
        show_path = k.show_pathway(str(kegg_path), keggid={str(kegg_id): "red"})
        res1 =  k.get(str(kegg_path)+"/image") # same as : res =  s.get("hsa05130","image")
        cf = open(str(uni_id)+"_path.png","wb")
        cf.write(res1)
    cf.close()
    print (cf, " saved!")

    bf.close()
    print (bf,"\n", " saved!")

    #map2 = k.mapping(fr='KEGG_ID', to='ACC', format='tab', query='hsa:7535')
    #for i in map2.values():
    #print i#'P43403'
    print ("\nBioservices - Kegg - ACCESS Gene - DONE!")
    #print f1, "\n",f2,"\n", f3,"\n", f4 ,"\n", " are saved!"
    theInputFiles.append(bf.name)
    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    print ("(complete!)")

#from_gene_id_get_kegg("LINJ_14_0700")
