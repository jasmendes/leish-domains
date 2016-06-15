# -*- coding: cp1252 -*-


# SUMMARY - Result Analysis, on 29-May-2015
import datetime, time
time_begin = datetime.datetime.fromtimestamp(time.time())


import sys, os
#print "FILE PATH:", sys.argv[0]
#print "FILE NAME:", os.path.basename(sys.argv[0])
#print sys.path
sys.path.append("C:\Python273\Lib\site-packages")
sys.path.append("C:\Python279\Lib\site-packages")
sys.path.append("C:\Python276\Lib\site-packages")
sys.path.append("C:\Python276\Lib")
#print sys.path

#from matplotlib import *
#from matplotlib_venn import *


import matplotlib_venn

from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles# module venn##""!EQW

import collections
from collections import Counter

import numbers
import decimal

from scipy import stats

from Bio import Entrez

import numpy as np
         

from scipy import *
import Tkinter as Tk

#from tk_tax_id_v154 import updatecombolist2
# when importing this it does not find functions on it. why?

from bioservices_functions_tk import *

#from auxiliary_functions_tk import *


import winsound

global theInputFiles


from pylab import figure, show, legend, ylabel
import winsound
########################################################################

########################################################################
##### MAKE IPR2GO DICT NAMES
#######################################################################3

kf = "_bmq_ipr_go.txt"

def get_dicts_ipr2go(kf):

    """Using the most common dict (COLLECTIONS) file from last functions,
    we are able to build a default DB set
    from  ipr2go EBI dictionary (kf)       
    with ipr:iprname, go:gonames  
                                  """


    #from os.path import basename
    # now you can call it directly with basename
    #print basename(af)
    #str(af[:-4])+
    global theInputFiles

    from collections import OrderedDict
    
    print "dict ipr2go = " , kf
    
    
    iprs1 = []
    iprs2 = []
    gos1 = []
    gos2 = []
    iprn_dict = {}
    gon_dict = {}

    FileCheck(kf)
    kfile = open(kf,"r")
    for i in kfile.readlines():
        i = i.strip()
        entry = i.split("\t")
        if len(entry) > 1:
            nameipr = entry[1]
            namego = entry[-2]
            for j in entry:
                if j[:3] == "IPR":
                    #print j
                    iprn_dict[j] = nameipr # [1]
                elif j[:3] == "GO:":
                    #print j
                    if len(namego) > 50:
                        namego = namego[:50]
                    gon_dict[j] = namego

    return iprn_dict, gon_dict



###############################################################################
#####     VENN2 AC vs AC
################################################################################


def read_tab_enrichment_acs(theFile1, theFile2):

    # Reads file as tab
    # Finds POSITION [0] AND COMPARES AvsB

    dict_iprn , dict_gon = get_dicts_ipr2go(kf)
    print "\nRequires the accessions to compare in position [0], both files ."

    
    af = open(theFile1, "r")
    bf = open(theFile2, "r")

    fields = []
    while len(fields) == 0:
        first = af.readline()
        first = first.strip()
        print first
        fields = first.split("\t")
        if len(fields) >1:
            break
 
    protdata = []

    iprs1 = []
    gos1 = []
    iprs2 = []
    gos2 = []
    uniq_iprs1 = []
    uniq_gos1 = []
    uniq_iprs2 = []
    uniq_gos2 = []

    ipr_names1 = []
    uniq_ipr_names1 = []

    ipr_names2 = []
    uniq_ipr_names2 = []

    acs1 = []
    acs2 = []
    for i in af.readlines():
        #print i
        i = i.strip()
        line = i.split("\t")
        acs1.append(line[0])
        for i in line:
            
            if "supermatch_entry_name" in fields: # IPR NAME
                pos = fields.index('supermatch_entry_name')
                #print pos
                ipr_name = line[pos] if len(line) > 2 else "null"
                ipr_names1.append(ipr_name)
                if ipr_name not in uniq_ipr_names1:
                    uniq_ipr_names1.append(ipr_name)
                    #print ipr_name
        for i in line:
            if i.startswith("IPR"): # IPR ID
                ipr1 = i
                #print i
                if len (ipr1) == 9:
                    #print ipr1, "oi"
                    iprs1.append(ipr1)
                    if ipr1 not in uniq_iprs1:
                        uniq_iprs1.append(ipr1)
                else:
                    #iprs1 = iprs1.split("; ")
                    for i in ipr1.split("; "):
                        #print i, " IPR FROM UNIPROT"
                        iprs1.append(i)
                        if i not in uniq_iprs1:
                            uniq_iprs1.append(i)
                    
        for i in line:
            if i.startswith("GO:"): # GO
                go = i
                #print go
                gos1.append(i)
                if i not in uniq_gos1:
                    uniq_gos1.append(i)

    for i in bf.readlines():
        #print i
        i = i.strip()
        line = i.split("\t")
        acs2.append(line[0])
        for i in line:
            if "supermatch_entry_name" in fields: # IPR NAME
                pos = fields.index('supermatch_entry_name')
                #print pos
                ipr_name = line[pos] if len(line) > 2 else "null"
                ipr_names2.append(ipr_name)
                if ipr_name not in uniq_ipr_names2:
                    uniq_ipr_names2.append(ipr_name)
                    #print ipr_name
        for i in line:
            if i.startswith("IPR"): # IPR ID
                ipr2 = i
                #print i
                if len (ipr2) == 9:
                    #print ipr1, "oi"
                    iprs2.append(ipr2)
                    if ipr2 not in uniq_iprs2:
                        uniq_iprs2.append(ipr2)
                else:
                    #iprs1 = iprs1.split("; ")
                    for i in ipr2.split("; "):
                        #print i, " IPR FROM UNIPROT"
                        iprs2.append(i)
                        if i not in uniq_iprs2:
                            uniq_iprs2.append(i)
                    
        for i in line:
            if i.startswith("GO:"): # GO
                go = i
                #print go
                gos2.append(go)
                if i not in uniq_gos2:
                    uniq_gos2.append(i)
    
     #f = ['a','b','d','c']
    one = set(acs1)
    two =set(acs2)
    xf = open(str(theFile1[:-4])+"_"+str(len(two))+"_venn2_acs.txt","w") #s = ['a','b','c']
    xf.write("VENN2  DIAGRAM \n"+str(theFile1)+"\n"+str(theFile2)+"\n")
   

    print "1 Intersection 2 = ", "\t", len(one.intersection(two))
    #print one.intersection(two)   

    

    #**set(['a', 'c', 'b'])**
    print "\nAll ACS ", af.name, bf.name
    print "1 Union 2 = ", "\t", len(one.union(two))
    ab_union = one.union(two)


   
    #**set(['a', 'c', 'b', 'd'])**
    print "\nACS ", af.name, bf.name
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - two.intersection(one))
    #**set(['d'])**

    #a = [1, 2, 3, 4, 5]
    #b = [9, 8, 7, 6, 5]
    #i for i, j in zip(a, b) if i == j]
    print "\nACS in common: ", af.name, bf.name
    print "2 Intersection 1 = ", "\t", len(two.intersection(one))
    #print one.intersection(two)

    #**set(['a', 'c', 'b'])**
    print "\nAll ACS ", af.name, bf.name
    print "2 Union 1 = ", "\t", len(two.union(one))
    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nACS ", af.name, bf.name
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(two.union(one)  - one.intersection(two))

    a = len(one)  - len(one.intersection(two))
    b = len(two)  - len(two.intersection(one))
    ab = len(one.intersection(two)) # len(one.union(two))
    
    aub = one.union(two)
    Ab_only = one - one.intersection(two)
    aB_only = two - two.intersection(one)
    ab_inter = one.intersection(two)

    

    xf.write("\n###################      Union AUB : "+ str(len(aub))+"\n\n")
    
    xf.write("\n######################      Only A : "+ str(len(Ab_only))+"\n")
    for i in Ab_only:
        if i in dict_iprn.keys():
            xf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            xf.write(str(i)+"\n")

    xf.write("\n######################      Only B : "+str(len(aB_only))+"\n")
    for i in aB_only:
        if i in dict_iprn.keys():
            xf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            xf.write(str(i)+"\n")
    xf.write("\n######################      Only AB : "+str(len(ab_inter))+"\n")
    for i in ab_inter:
        if i in dict_iprn.keys():
            xf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            xf.write(str(i)+"\n")
    #get_venn22(a,b,ab)

    #theInputFiles.append(xf.name)
    theInputFiles.append(xf.name)

    print a, b, ab

    ####        VENN  DIAGRAM         ######

    fig = plt.figure()
    fig.canvas.set_window_title(str(xf.name[:-4]))

    
    # Subset sizes
    s = (
        a,  # Ab
        b,  # aB
        ab,  # AB
        )

    v = venn2(subsets=s, set_labels=('A', 'B'))

    # Subset labels
    v.get_label_by_id('10').set_text(str(a))
    v.get_label_by_id('01').set_text(str(b))
    v.get_label_by_id('11').set_text(str(ab))

    # Subset colors
    v.get_patch_by_id('10').set_color('c')
    v.get_patch_by_id('01').set_color('#993333')
    v.get_patch_by_id('11').set_color('blue')

    # Subset alphas
    v.get_patch_by_id('10').set_alpha(0.4)
    v.get_patch_by_id('01').set_alpha(0.5)
    v.get_patch_by_id('11').set_alpha(0.7)

    
    plt.title(str(theFile1)+" vs.\n "+str(theFile2))

    # Border styles
    c = venn2_circles(subsets=s, linestyle='solid')
    c[0].set_ls('dashed')  # Line style
    c[0].set_lw(2.0)       # Line width

    
    plt.show()

    

    print "(complete!)"


#read_tab_enrichment_ipr("output/venns/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt",
                         #"output/venns/9606_ac_entry_ipr_name_gene_go_path_local_exi.txt")
#read_tab_enrichment_ipr("output/2pep_246_down_filter_246_tax5671_bs_uniprot.txt" ,"output/2pep_260_filter_260_tax5671_bs_uniprot.txt")
#read_tab_enrichment_ipr("output/5671_bmq_ac_ipr_go.txt", "output/5661_bmq_ac_ipr_go.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")
###read_tab_enrichment_ipr("output/tax5671_bmq_ac_ipr_gos.txt", "output/tax5671_bmq_ac_ipr.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")
#read_tab_enrichment_ipr("output/2peps_alls_filter_3998_tax5671_bs_uniprot.txt", "vy.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")
###############################################################################

#############################################################################################################################################################3
##################################         VENN 2    VENN 3   ################3
#######################################################3#################3


#import summary_plot_tk
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn2_circles


######################################     VENN 2     ####################################################

def read_tab_enrichment_ipr(theFile1, theFile2):

    # Reads file as tab
    # Finds IPR
    # Finds True and False s

    dict_iprn , dict_gon = get_dicts_ipr2go(kf)

    
    af = open(theFile1, "r")
    bf = open(theFile2, "r")

    fields = []
    while len(fields) == 0:
        first = af.readline()
        first = first.strip()
        print first
        fields = first.split("\t")
        if len(fields) >1:
            break
 
    protdata = []

    iprs1 = []
    gos1 = []
    iprs2 = []
    gos2 = []
    uniq_iprs1 = []
    uniq_gos1 = []
    uniq_iprs2 = []
    uniq_gos2 = []

    ipr_names1 = []
    uniq_ipr_names1 = []

    ipr_names2 = []
    uniq_ipr_names2 = []

    for i in af.readlines():
        #print i
        i = i.strip()
        line = i.split("\t")
        for i in line:
            if "supermatch_entry_name" in fields: # IPR NAME
                pos = fields.index('supermatch_entry_name')
                #print pos
                ipr_name = line[pos] if len(line) > 2 else "null"
                ipr_names1.append(ipr_name)
                if ipr_name not in uniq_ipr_names1:
                    uniq_ipr_names1.append(ipr_name)
                    #print ipr_name
        for i in line:
            if i.startswith("IPR"): # IPR ID
                ipr1 = i
                #print i
                if len (ipr1) == 9:
                    #print ipr1, "oi"
                    iprs1.append(ipr1)
                    if ipr1 not in uniq_iprs1:
                        uniq_iprs1.append(ipr1)
                else:
                    #iprs1 = iprs1.split("; ")
                    for i in ipr1.split("; "):
                        #print i, " IPR FROM UNIPROT"
                        iprs1.append(i)
                        if i not in uniq_iprs1:
                            uniq_iprs1.append(i)
                    
        for i in line:
            if i.startswith("GO:"): # GO
                go = i
                #print go
                gos1.append(i)
                if i not in uniq_gos1:
                    uniq_gos1.append(i)

    for i in bf.readlines():
        #print i
        i = i.strip()
        line = i.split("\t")
        for i in line:
            if "supermatch_entry_name" in fields: # IPR NAME
                pos = fields.index('supermatch_entry_name')
                #print pos
                ipr_name = line[pos] if len(line) > 2 else "null"
                ipr_names2.append(ipr_name)
                if ipr_name not in uniq_ipr_names2:
                    uniq_ipr_names2.append(ipr_name)
                    #print ipr_name
        for i in line:
            if i.startswith("IPR"): # IPR ID
                ipr2 = i
                #print i
                if len (ipr2) == 9:
                    #print ipr1, "oi"
                    iprs2.append(ipr2)
                    if ipr2 not in uniq_iprs2:
                        uniq_iprs2.append(ipr2)
                else:
                    #iprs1 = iprs1.split("; ")
                    for i in ipr2.split("; "):
                        #print i, " IPR FROM UNIPROT"
                        iprs2.append(i)
                        if i not in uniq_iprs2:
                            uniq_iprs2.append(i)
                    
        for i in line:
            if i.startswith("GO:"): # GO
                go = i
                #print go
                gos2.append(go)
                if i not in uniq_gos2:
                    uniq_gos2.append(i)
    
     #f = ['a','b','d','c']
    one = set(uniq_iprs1)
    two =set(uniq_iprs2)
    xf = open(str(theFile1[:-4])+"_"+str(len(two))+"_venn2.txt","w") #s = ['a','b','c']
    xf.write("VENN2  DIAGRAM \n"+str(theFile1)+"\n"+str(theFile2)+"\n")
   

    print "UNIQUE IPRS 1 vs IPRS2 ", len(uniq_iprs1), len(uniq_iprs2)
    print "TOTAL IPRS 1 vs IPRS2 ", len(iprs1), len(iprs2)
    print "\nDomains in common: ", af.name, bf.name
    print "1 Intersection 2 = ", "\t", len(one.intersection(two))
    #print one.intersection(two)   

    

    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", af.name, bf.name
    print "1 Union 2 = ", "\t", len(one.union(two))
    ab_union = one.union(two)


   
    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", af.name, bf.name
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - two.intersection(one))
    #**set(['d'])**

    #a = [1, 2, 3, 4, 5]
    #b = [9, 8, 7, 6, 5]
    #i for i, j in zip(a, b) if i == j]
    print "\nDomains in common: ", af.name, bf.name
    print "2 Intersection 1 = ", "\t", len(two.intersection(one))
    #print one.intersection(two)

    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", af.name, bf.name
    print "2 Union 1 = ", "\t", len(two.union(one))
    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", af.name, bf.name
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(two.union(one)  - one.intersection(two))

    a = len(one)  - len(one.intersection(two))
    b = len(two)  - len(two.intersection(one))
    ab = len(one.intersection(two)) # len(one.union(two))
    
    aub = one.union(two)
    Ab_only = one - one.intersection(two)
    aB_only = two - two.intersection(one)
    ab_inter = one.intersection(two)

    

    xf.write("\n###################      Union AUB : "+ str(len(aub))+"\n\n")
    
    xf.write("\n######################      Only A : "+ str(len(Ab_only))+"\n")
    for i in Ab_only:
        if i in dict_iprn.keys():
            xf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            xf.write(str(i)+"\n")

    xf.write("\n######################      Only B : "+str(len(aB_only))+"\n")
    for i in aB_only:
        if i in dict_iprn.keys():
            xf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            xf.write(str(i)+"\n")
    xf.write("\n######################      Only AB : "+str(len(ab_inter))+"\n")
    for i in ab_inter:
        if i in dict_iprn.keys():
            xf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            xf.write(str(i)+"\n")
    #get_venn22(a,b,ab)

    #theInputFiles.append(xf.name)
    theInputFiles.append(xf.name)

    print a, b, ab

    ####        VENN  DIAGRAM         ######

    fig = plt.figure()
    fig.canvas.set_window_title(str(xf.name[:-4]))

    
    # Subset sizes
    s = (
        a,  # Ab
        b,  # aB
        ab,  # AB
        )

    v = venn2(subsets=s, set_labels=('A', 'B'))

    # Subset labels
    v.get_label_by_id('10').set_text(str(a))
    v.get_label_by_id('01').set_text(str(b))
    v.get_label_by_id('11').set_text(str(ab))

    # Subset colors
    v.get_patch_by_id('10').set_color('c')
    v.get_patch_by_id('01').set_color('#993333')
    v.get_patch_by_id('11').set_color('blue')

    # Subset alphas
    v.get_patch_by_id('10').set_alpha(0.4)
    v.get_patch_by_id('01').set_alpha(0.5)
    v.get_patch_by_id('11').set_alpha(0.7)

    
    plt.title(str(theFile1)+" vs.\n "+str(theFile2))

    # Border styles
    c = venn2_circles(subsets=s, linestyle='solid')
    c[0].set_ls('dashed')  # Line style
    c[0].set_lw(2.0)       # Line width

    
    plt.show()

    

    print "(complete!)"


#read_tab_enrichment_ipr("output/venns/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt",
                         #"output/venns/9606_ac_entry_ipr_name_gene_go_path_local_exi.txt")
#read_tab_enrichment_ipr("output/2pep_246_down_filter_246_tax5671_bs_uniprot.txt" ,"output/2pep_260_filter_260_tax5671_bs_uniprot.txt")
#read_tab_enrichment_ipr("output/5671_bmq_ac_ipr_go.txt", "output/5661_bmq_ac_ipr_go.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")
###read_tab_enrichment_ipr("output/tax5671_bmq_ac_ipr_gos.txt", "output/tax5671_bmq_ac_ipr.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")
#read_tab_enrichment_ipr("output/2peps_alls_filter_3998_tax5671_bs_uniprot.txt", "vy.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")



    

########################################################
## VENN 3
############################################################

import pylab as plt
from matplotlib_venn import venn3, venn3_circles



def get_domain_venn_diagram1(theFilename1, theFilename2, theFilename3): #, theFilename4, theOutput):

    "Entry  Entry name  Protein names   InterPro    Gene ontology IDs   Gene ontology (GO)  Gene names  Pathway Protein existence   Organism    Organism ID"
    "protein_accession	protein_name	entry_ac	entry_name	entry type	go_id	go_term	go_root	protein_database_name	tax_name	tax_id"

    global theInputFiles

    dict_iprn , dict_gon = get_dicts_ipr2go(kf)
    
    # Given this format will produce countings for all InterPro entries and GO terms
    aFile = open(theFilename1, "rU")
    bFile = open(theFilename2, "rU")
    cFile = open(theFilename3, "rU")
    #dFile = open(theFilename4, "rU") # Human
    #eFile = open(theOutput, "wt")
    Lin_domains = {} # 5148
    Ldo_domains = {} # 2893
    Lma_domains = {} # 5276
    Hum_domains= {} # 11748
    Lin_unique_domains = []
    Lma_unique_domains = []
    Ldo_unique_domains = []
    Hum_unique_domains = []
    count = 0
    count1 = 0
    count2 = 0
    count_not = 0
    count_in_ma = 0
    count_in_do = 0 
    first1 = aFile.readline()
    fields1 = first1.split("\t")
    first2 = bFile.readline()
    fields2 = first1.split("\t")
    first3 = cFile.readline()
    fields3 = first1.split("\t")
    for aRow in aFile.readlines():
        #print aRow
        aProtein_entry = aRow.strip().split("\t")
        if len(aProtein_entry) == 0:
            continue

        elif len(aProtein_entry[3])== 0:
            continue
        elif len(aProtein_entry[3])== 9:
            Lin_domains[aProtein_entry[3]] = aRow
            if aProtein_entry[3] not in Lin_unique_domains:
                Lin_unique_domains.append(aProtein_entry[3])
        else:
            for j in aProtein_entry:
                if j[:3] == "IPR":
                    iprs3 = j
                    for i in iprs3.split("; "):
                        Lin_domains[i] = aRow
                        if i not in Lin_unique_domains:
                            Lin_unique_domains.append(i)
            
    
    for aRow in bFile.readlines():
        #print aRow
        aProtein_entry = aRow.strip().split("\t")
        if len(aProtein_entry) == 0:
            continue
        elif len(aProtein_entry[3])== 0:
            continue
        elif len(aProtein_entry[3])== 9:
            Lma_domains[aProtein_entry[3]] = aRow
            if aProtein_entry[3] not in Lma_unique_domains:
                Lma_unique_domains.append(aProtein_entry[3])
        else:
            for j in aProtein_entry:
                if j[:3] == "IPR":
                    iprs2 = j
                    for i in iprs2.split("; "):
                        Lma_domains[i] = aRow
                        if i not in Lma_unique_domains:
                            Lma_unique_domains.append(i)
    
    for aRow in cFile.readlines():
        #print aRow
        aProtein_entry = aRow.strip().split("\t")
        if len(aProtein_entry) == 0:
            continue
        elif len(aProtein_entry[3])== 0:
            continue
        elif len(aProtein_entry[3])== 9:
            Ldo_domains[aProtein_entry[3]] = aRow
            if aProtein_entry[3] not in Ldo_unique_domains:
                Ldo_unique_domains.append(aProtein_entry[3])
        else:
            for j in aProtein_entry:
                if j[:3] == "IPR":
                    iprs3 = j
                    for i in iprs3.split("; "):
                        #print i
                        Ldo_domains[i] = aRow
                        if i not in Ldo_unique_domains:
                            Ldo_unique_domains.append(i)
                
    #for i in Lin_domains and Lma_domains:
        #print "Lin:", str(Lin_domains[i])
        #print "Lma:", str(Lma_domains[i])
     #   count_in_ma = count_in_ma + 1
        #3elif i in Lma_domains[i] and not Ldo_domains[i]:
            #count_in_ma = count_in_ma + 1
            #print "Lin U Lma:", i
        #elif i in Ldo_domains[i] and not Lma_domains[i]:
            #count_in_do = count_in_do + 1
        #elif i not in Lma_domains and not Ldo_domains:
            #print "NOT", "Linfantum:", Lin_domains[i]
            #count_not = count_not + 1
        #else:
            #print "WHY"
    #print "Linfantum domains in Lma and Ldo", count #(5243 Lma) 4718 Lma in Lin , 4827 Lin in Lma

    print len(Lin_domains), "\t in ", aFile.name # 1561
    print len(Lma_domains), "\t in ", bFile.name  # 1570
    print len(Ldo_domains), "\t in ", cFile.name  # 1069
    #print len(Lin_unique_domains), "\t in ", aFile.name  #1561
    #print len(Lma_unique_domains), "\t in ", bFile.name #1570
    #print len(Ldo_unique_domains), "\t in ", cFile.name  #1069
    
    """for aRow in dFile.readlines():
        #print aRow
        aProtein_entry = aRow.strip().split("\t")
        if len(aProtein_entry[3])== 0:
            continue
        elif len(aProtein_entry[3])== 9:
            Hum_domains[aProtein_entry[3]] = aRow
            if aProtein_entry[3] not in Hum_unique_domains:
                Hum_unique_domains.append(aProtein_entry[3])
        else:
            for j in aProtein_entry:
                if j[:3] == "IPR":
                    iprs3 = j
                    for i in iprs3.split("; "):
                        #print i
                        Hum_domains[i] = aRow
                        if i not in Hum_unique_domains:
                            Hum_unique_domains.append(i)
    print len(Hum_unique_domains), "\t in ", dFile.name  #1069"""



    count_3 = 0
    count_o1 = 0
    count_o2 = 0
    count_n1 = 0
    count_n2 = 0
    count_no = 0
    count_ma = 0
    count_do = 0
    
    for i in Lin_domains.keys(): # 3610
        if i in Ldo_domains.keys() and Lma_domains.keys():
            count_3 += 1 # 3548
        if i in Ldo_domains.keys() and not Lma_domains.keys(): # 3607
            count_o1 += 1 # 0
        if i in Ldo_domains.keys():
            count_do += 1
        if i in Lma_domains.keys() and not Ldo_domains.keys(): # 3572
            count_o2 += 1 # 0
        if i in Lma_domains.keys():
            count_ma += 1
        if i not in Ldo_domains.keys(): 
            count_n1 += 1 # 48
            #print i, Lin_domains[i]
        if i not in Lma_domains.keys():
            count_n2 += 1 # 62
            #print i, Lin_domains[i]
        if i not in Lma_domains.keys() and not Ldo_domains.keys():
            count_no += 1 # 0
            #print i, Lin_domains[i]

    print count_3 # 2720 total lin
    print count_o1 # 55 lin exclusive
    print count_o2
    print count_n1, count_n2, count_no# 0 in Ldo, 74 in Lma
    print count_do, count_ma # 3562, 3548
    """3548
0
0
62 48 0"""
    
    one = set(Lin_unique_domains)
    two =set(Lma_unique_domains)
    three = set(Ldo_unique_domains)
    #four = set(Hum_unique_domains)
    
    print "\nDomains in common: ", aFile.name, bFile.name
    print "1 Intersection 2 = ", "\t", len(one.intersection(two))
    #print one.intersection(two)

    print "\nAll Domains ", aFile.name, bFile.name
    print "1 Union 2 = ", "\t", len(one.union(two))
    
    print "\nDomains ", aFile.name, bFile.name
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - one.intersection(two))

    #######
    
    print "\nDomains in common: ", bFile.name, cFile.name
    print "2 Intersection 1 = ", "\t", len(two.intersection(three))
    #print one.intersection(two)

    print "\nAll Domains ", bFile.name, cFile.name
    print "2 Union 1 = ", "\t", len(two.union(three))
    
    print "\nDomains ", bFile.name, cFile.name
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(two.union(three)  - two.intersection(three))

    #######
    print "\nDomains in common: ", aFile.name, cFile.name
    print "2 Intersection 1 = ", "\t", len(one.intersection(three))
    este1 = one.intersection(three)
    este2 = one.intersection(three)
    este3 = one.intersection(three)
    
    #print one.intersection(two)

    print "\nAll Domains ", aFile.name, cFile.name
    print "2 Union 1 = ", "\t", len(one.union(three))
    
    print "\nDomains ", aFile.name, cFile.name
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(one.union(three)  - one.intersection(three))

    #######
    print "\nDomains in common: ", aFile.name, cFile.name
    print "2 Intersection 1 = ", "\t", len(one.intersection(three))
    #print one.intersection(two)

    print "\nAll Domains ", aFile.name, cFile.name
    print "2 Union 1 = ", "\t", len(one.union(three))
    
    print "\nDomains ", aFile.name, cFile.name
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(one.union(three)  - one.intersection(three))



    
    """####### Human ####
    print "\nDomains in common: ", aFile.name, dFile.name
    print "2 Intersection 1 = ", "\t", len(one.intersection(four))
    
    look = list(one.intersection(two))
    print look[0:25]

    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", aFile.name, dFile.name
    print "2 Union 1 = ", "\t", len(one.union(four))
    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", aFile.name, dFile.name
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(one.union(four)  - one.intersection(four))"""

    va = one - one.intersection(two.union(three)) # only A
    vb = two  - two.intersection(one.union(three)) # only B
    vc = three - three.intersection(one.union(two)) # Only C
    

    abc_inter = set.intersection(one,two,three) # only ABC
    
    abc_union = set.union(one,two,three)

    ab_inter = set.intersection(one,two)
    bc_inter = set.intersection(two,three)
    ac_inter = set.intersection(two,three)

    vbc = two.intersection(three) - two.intersection(three.intersection(one)) #ab_inter - abc_inter  # only AB
    vab = one.intersection(two) - one.intersection(two.intersection(three)) #bc_inter - abc_inter  # only BC
    vac = one.intersection(three) - one.intersection(three.intersection(two))

    
    Abc = len(one) - len(one.intersection(two.union(three)))
    aBc = len(two)  - len(two.intersection(one.union(three)))
    abC = len(three) - len(three.intersection(one.union(two)))
    #ab = len(one.intersection(two)) # len(one.union(two))

    aBC = len(two.intersection(three)) - len(two.intersection(three.intersection(one)))
    ABc = len(one.intersection(two)) - len(one.intersection(two.intersection(three)))
    AbC = len(one.intersection(three)) - len(one.intersection(three.intersection(two)))
    ABC = len(set.intersection(one,two,three))

    
    
    
    soma = Abc + aBc + abC + ABc + AbC + aBC + ABC
    print Abc ,aBc, abC, ABc,AbC,aBC, ABC
    print Abc, aBc, abC, len(va), len(vb), len(vc)
    
    print "Sum: ", soma
    print "\n TOTAL UNIQ IPRS 1 vs 2 vs 3"
    print len(one), len(two), len(three)

    print "Intersection ABC :", len(abc_inter)
    print "Union ABC : ", len(abc_union)


    
    # Create Title Name

    ven_names = ()

    theinputs= []
    theinputs.append(theFilename1)
    theinputs.append(theFilename2)
    theinputs.append(theFilename3)

    for i in theinputs:

        if "/" in i:
            i = i.split("/")
            name1 = i[1]
            ven_name = name1.split("_")[0]
            ven_names = ven_names + (ven_name,)
        elif "\\" in i:
            i = i.split("\\")
            name1 = i[1]
            ven_name = name1.split("_")[0]
            ven_names = ven_names + (ven_name,)
        else:
            bar_name = name1.split("_")[0]
            ven_names = ven_names + (ven_name,)


      
    yf = open(theFilename1[:-4]+"_"+ven_names[1]+"_"+ven_names[2]+"_"+str(ABC)+"_venn3.txt","w")
    
    yf.write("VENN3  DIAGRAM \n"+str(theFilename1)+"\n"+str(theFilename2)+"\n"+str(theFilename3)+"\n")
    
    yf.write("\n###################      Union AUBUC : "+ str(len(abc_union))+"\n\n")
    yf.write("\n###################     Only A : "+ str(Abc)+"\n")
    for i in va:
        if i in dict_iprn.keys():
            yf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            yf.write(str(i)+"\n")

    yf.write("\n###################     Only B : "+str(aBc)+"\n")
    for i in vb:
        if i in dict_iprn.keys():
            yf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            yf.write(str(i)+"\n")

    yf.write("\n###################     Only C : "+str(abC)+"\n")
    for i in vc:
        if i in dict_iprn.keys():
            yf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            yf.write(str(i)+"\n")
    yf.write("\n###################     Only AC : "+str(AbC)+"\n")
    for i in vac:
        if i in dict_iprn.keys():
            yf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            yf.write(str(i)+"\n")

    yf.write("\n###################     Only BC : "+str(aBC)+"\n")
    for i in vbc:
        if i in dict_iprn.keys():
            yf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            yf.write(str(i)+"\n")

    yf.write("\n###################     Only AB : "+str(ABc)+"\n")
    for i in vab:
        if i in dict_iprn.keys():
            yf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            yf.write(str(i)+"\n")

    yf.write("\n###################     Intersection  ABC : "+str(ABC)+"\n")
    for i in abc_inter:
        if i in dict_iprn.keys():
            yf.write(str(i)+"\t"+str(dict_iprn[i])+"\n")
        else:
            yf.write(str(i)+"\n")
        
    #get_venn22(a,b,ab)
    #get_venn22(a,b,ab)
    

    #abc_union = set.union(one,two,three)
    
    

    
    ###################   VENN3 DIAGRAM    #########################
    # Subset sizes
    s = (
        Abc,    # Abc
        aBc,    # aBc
        ABc,    # ABc
        abC,    # abC
        AbC,    # AbC
        aBC,  # aBC
        ABC,    # ABC
        )


    ##### VENN 3  ###

    fig = plt.figure()
    fig.canvas.set_window_title(str(yf.name[:-4]))


    v = venn3(subsets=s, set_labels=ven_names)

    # Subset labels
    v.get_label_by_id('100').set_text(str(Abc))
    v.get_label_by_id('010').set_text(str(aBc))
    v.get_label_by_id('110').set_text(str(ABc))
    v.get_label_by_id('001').set_text(str(abC))
    v.get_label_by_id('101').set_text(str(AbC))
    v.get_label_by_id('011').set_text(str(aBC))
    v.get_label_by_id('111').set_text(str(ABC))

    # Subset colors
    v.get_patch_by_id('100').set_color('c')
    v.get_patch_by_id('010').set_color('#993333')
    v.get_patch_by_id('110').set_color('blue')

    # Subset alphas
    v.get_patch_by_id('101').set_alpha(0.4)
    v.get_patch_by_id('011').set_alpha(1.0)
    v.get_patch_by_id('111').set_alpha(0.7)

    # Border styles
    c = venn3_circles(subsets=s, linestyle='solid')
    c[0].set_ls('dotted')  # Line style
    c[1].set_ls('dashed')
    c[2].set_lw(1.0)       # Line width

    plt.title(str(theFilename1)+" vs.\n "+str(theFilename2)+" vs.\n "+str(theFilename3))

    
    plt.show()

    yf.close()
    print yf.name, " Saved!"
    

    theInputFiles.append(yf.name)

    print "(complete!)"


#get_domain_venn_diagram1("output/venns/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt",
           #             "output/venns/5664_ac_entry_ipr_name_gene_go_path_local_exi.txt",
            #           "output/venns/5661_ac_entry_ipr_name_gene_go_path_local_exi.txt")
             #        #    "output/venns/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt", 
             #           "output/venns/Lin_complex_venn.txt") # "output/venns/9606_ac_entry_ipr_name_gene_go_path_local_exi.txt"
#get_domain_venn_diagram1("20150917_VENN2_4FC/all4FCup_filter_733_all_FC_filter_5213_tax5671_bs_uniprot.txt","20150917_VENN2_4FC/all202FC_filter_4161_all_FC_filter_5213_tax5671_bs_uniprot.txt","20150917_VENN2_4FC/all4FCdown_filter_322_all_FC_filter_5213_tax5671_bs_uniprot.txt")

#################################################################################
#####   SCATTER
################################################################################

import numpy as np
import matplotlib.pyplot as plt

# http://alstatr.blogspot.pt/2013/11/python-venn-diagram.html

def scatter_legend(sf):

    N = 10
    data = np.random.random((N, 4))
    labels = ['point{0}'.format(i) for i in range(N)]
    plt.subplots_adjust(bottom = 0.1)
    plt.scatter(
        data[:, 0], data[:, 1], marker = 'o', c = data[:, 2], s = data[:, 3]*1500,
        cmap = plt.get_cmap('Spectral'))

    for label, x, y in zip(labels, data[:, 0], data[:, 1]):
        plt.annotate(
            label,
            xy = (x, y), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

    plt.show()
    
#scatter_legend()


##############################################################################3

####################################################################
#### PIE
##################################################################
    


def count_go_pie(af, kf): # n = 1000

    """Using the most common dict (COLLECTIONS) file from last functions,
    we are able to build a default proteome set (bf) 
    and compare it with input proteome set (af) 
    combining de ipr2go EBI dictionary (kf)                                     """


    #from os.path import basename
    # now you can call it directly with basename
    #print basename(af)
    #str(af[:-4])+
    global theInputFiles

    from collections import OrderedDict
    
    print "choice 1 = ", af
    
    print "dict ipr2go = " , kf
    
    cfile = str(af[:-4])+"_counts_ipr_go.txt"
    cf = open(cfile, "w")
    
    iprs1 = []
    iprs2 = []
    gos1 = []
    gos2 = []
    iprn_dict = {}
    gon_dict = {}

    FileCheck(kf)
    kfile = open(kf,"r")
    fields = kfile.readline()
    for i in kfile.readlines():
        i = i.strip()
        entry = i.split("\t")
        if len(entry) > 1:
            nameipr = entry[1]
            namego = entry[-2]
            for j in entry:
                if j[:3] == "IPR":
                    #print j
                    iprn_dict[j] = nameipr # [1]
                elif j[:3] == "GO:":
                    #print j
                    if len(namego) > 50:
                        namego = namego[:50]
                    gon_dict[j] = namego
                    


    for i in kfile.readlines():
        i = i.strip()
        entry = i.split("\t")
        #print entry
        if len(entry) > 3 :
            namego = entry[-2]
            #print namego
            for j in entry:
                if j[:3] == "GO:":
                    #print j
                    gon_dict[j] = namego
    
    acs1 = []

    iacs1 = []
    iprc1 = 0
    goc1 = 0
    iprs11 = []
    gos11 = []
    all_ac1 = []
    all_ac2 = []
    gacs1 = []
    gacs2 = []
    notfound = 0
    num_lines1 = sum(1 for line in open(af))
    with open(af,"r") as f1:
        for i in f1.readlines()[1:]:
            #print i
            i = i.strip()
            items = i.split("\t")
            ac = items[0]
            if ac not in all_ac1:
                all_ac1.append(ac)
            for j in items:
                #print j
                if j[0:3] == "IPR":
                    iprc1 += 1
                    if ac not in iacs1:
                        iacs1.append(ac)                    
                    #print index(j)
                    #print j.split("; ")
                    for i in j.split("; "):
                        iprs1.append(i)
                        if i not in iprs11:
                            iprs11.append(i)
                elif j[0:3] == "GO:":
                    goc1 += 1
                    #print j.split("; ")
                    if ac not in gacs1:
                        gacs1.append(ac)
                    for i in j.split("; "):
                        gos1.append(i)
                        if i not in gos11:
                            gos11.append(i)
                #else:
                #notfound += 1
                
                

    iacs2 = []
    #acs2 = []
    iprc2 = 0
    goc2 = 0
    iprs22 = []
    gos22 = []

    no_go = len(all_ac1) - len(gacs1) 
    no_ipr = len(all_ac1) - len(iacs1) 

    #print no_go, no_ipr
    #print len(gacs1), len(iacs1)
    #print 
    
                

    

    ##############3   PIE CHART     #######################3
    
    
    
    c  = Counter(iprs1)
    print "\nTOTAL iprs1 found : ", len(iprs1) # 6597
    print "UNIQUE iprs1 : ", len(c), "\n", c.most_common(10) # 2088

    cf.write("########################## UNIQUE IPRS: "+str(len(c))+"\n")

    c_sorted_by_value = OrderedDict(sorted(c.items(), key=lambda x: x[1]), reverse = True)
    for i in c_sorted_by_value:
        if i[:3] == "IPR":
            if i in iprn_dict.keys():
                cf.write(str(i)+"\t"+str(c[i])+"\t"+str(iprn_dict[i])+"\n")
            else:
                cf.write(str(i)+"\t"+str(c[i])+"\n")            

            
            # print i, f[i]
        
    #for i in dict(c.most_common()):
     #   cf.write(str(i)+"\t"+str(c[i])+"\n")
        # print i, c[i]
        
    ### COVERAGE IPR 1
    #print iprc1, acs1
    try:
        covi =  (len(iacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enri = (len(iprs11)/float(len(iprs1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation IPR 1 :\n\tAc Coverage: ", ("{0:.2%}".format(covi))#
        print "\tEnrichment", ("{0:.2%}".format(enri))
    except:
        print "No IPRs 1 found"



    ### GO SUM ###
    ### COVERAGE GO 1

   
    
    f  = Counter(gos1)
    print "\nTOTAL GO IDS 1 found : ", len(gos1) # 6597
    print "UNIQUE GO IDs 1: ", len(f), "\n", f.most_common(10) # 1183

    cf.write("\n\n######################## UNIQUE GOS: "+str(len(f))+"\n")

    f_sorted_by_value = OrderedDict(sorted(f.items(), key=lambda x: x[1]), reverse=True)
    for i in f_sorted_by_value:
        if i[:3] == "GO:":
            if i in gon_dict.keys():
                cf.write(str(i)+"\t"+str(f[i])+"\t"+str(gon_dict[i])+"\n")
            else:
                cf.write(str(i)+"\t"+str(f[i])+"\n")
                # print i, f[i]

    try:
        covg = (len(gacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enrg = (len(gos11)/float(len(gos1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 1 :\n\tAc Coverage : ", ("{0:.2%}".format(covg))#
        print "\tEnrichment", ("{0:.2%}".format(enrg))
    except:
        print "No GOs 1 found!"

        
    #pie_plot_legends(ipr_pie_names, ipr_pie_nums)
    #pie_plot_legends(go_pie_names, go_pie_nums)
    # The slices will be ordered and plotted counter-clockwise.

    theInputFiles.append(cf.name)

    go_pie_names = []
    go_pie_names2 = []
    go_pie_nums = []
    
    fout = 0
    repercentgo = 0

    for i, j in f.most_common():
        if i[:3] == "GO:":
            go_pie_names.append(i)
            go_pie_nums.append(j)
            new = str(j) 
            if i in gon_dict.keys():
                new = str(j) +" "+ gon_dict[i]
            
            go_pie_names2.append(new) # TO LABELS
            fout = fout + j 
            #print fout, len(iprs1), repercent
            repercentgo = fout / len(gos1)
            
            #half = int(len(gos1)*0.3)
            half = int(len(gos1)-fout)
            if len(go_pie_names) > 13 :
                go_pie_names.append("Other")
                go_pie_names.append("No GO ID")
                go_pie_nums.append(half)
                
                new = str(half) +" Other"
                
                go_pie_names2.append(new) # TO LABELS
                
                if no_go<0:
                    go_pie_nums.append(no_go)
                    new2 = str(no_go) + " No GO ID"
                    go_pie_names2.append(new2) # TO LABELS
                break


    repercent2go = fout / (len(gos11)+1)
    print repercentgo, repercent2go
        

    

    ipr_pie_names = []
    ipr_pie_nums = []

    ipr_pie_names2 = []

    fout = 0
    repercentipr = 0
    for i, j in c.most_common():
        if i[:3] == "IPR":
            ipr_pie_names.append(i)
            ipr_pie_nums.append(j)
            new = str(j) +" "+ iprn_dict[i]
            ipr_pie_names2.append(new) # TO LABELS
            fout = fout + j 
            #print fout, len(iprs1), repercent
            repercentipr = fout / len(iprs11)
            if repercentipr > 0.5 :
                break
            
    repercent2ipr = fout / (len(iprs11)+1)
    #print repercentipr, repercent2ipr

    fig = plt.figure()
    fig.canvas.set_window_title(str(cf.name[:-4])) # no .name

    
    ## IPR ##
    #labels = ipr_pie_names # 'Frogs', 'Hogs', 'Dogs', 'Logs'
    #labels2 = ipr_pie_names2
    #sizes = ipr_pie_nums # [15, 30, 45, 10]

    ## GO ##
    labels = go_pie_names
    labels2 = go_pie_names2
    sizes = go_pie_nums 
    
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    explode = (2, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90)#, labels=labels
    # Set aspect ratio to be equal so that pie is drawn as a circle.

    plt.axis('equal')

    #plt.title(str(af))
    
    plt.legend(labels2, loc='best', shadow=True)
    #pietitle = "ACsGO=" +str(len(gacs1))+", TotGOs="+str(len(gos1))
    pietitle = " ACsGO= %s , TotGOs= %s , UniqGOs = %s"  % (len(gacs1),len(gos1), len(gos11))
    pietitle2 = "\nACs= %s, Coverage = %s , Enrichment = %s " % (len(all_ac1), "{0:.2%}".format(covg), "{0:.2%}".format(enrg))

    
    plt.suptitle("20 Most Common GO IDs in Protein Accessions " + pietitle+ pietitle2)
    #plt.legend(title="técnica")

    #ax.legend([extra, bar_0_10, bar_10_100], ("My explanatory text", "0-10", "10-100"))

    #plt.savefig(name)
    #print "Representing %s of the total IPRs " % repercentipr

    repercentgo = repercentgo * 100
    centgo = "{0:.1f}".format(repercentgo)

    repercent2ipr= repercent2ipr * 100
    centipr = "{0:.1f}".format(repercent2ipr)

    
    
    #print "Representing %s %% of the total IPRs " % centipr
    #print "Representing %s %% of the total GOs " % centgo

   

    plt.show()

    
    

    print "(complete!)"

    """## GO ##

    labels2 = go_pie_names # 'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes2 = go_pie_nums # [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    #explode = (0, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes2, labels=labels2, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90)
    # Set aspect ratio to be equal so that pie is drawn as a circle.

    plt.axis('equal')

    plt.show()"""

    print "\nAC with IPR : ", len(iacs1)
    print "TOTAL IPRs1 : ",len(iprs1)
    print "UNIQUE IPRs1 : ",len(iprs11)
    

    print "AC with GO : ", len(gacs1)
    print "TOTAL GOs1 : ",len(gos1)
    print "UNIQUE GOs1 : ",len(gos11)


    print"\nTOTAL UNIQUE ACS 1 :",len(all_ac1)
    print"\nTOTAL UNIQUE ACS 1 :",len(all_ac2)

        

    print "\n", cf.name, " saved!"

    
    
    #print ef, " saved!"
    #ef.close()
    #cff.close()
    #dff.close()
    kfile.close()
    cf.close()


#count_go_pie("250_down_all_filter_250_tax5671_bmq_ac_ipr_go.txt", "output/_bmq_ipr_go.txt")

#############################     IPR   PIE      ##########################

def count_ipr_pie(af, kf): # n = 1000

    """Using the most common dict (COLLECTIONS) file from last functions,
    we are able to analyze input proteome set (af) 
    combining de ipr2go EBI dictionary (kf)                                     """


    #from os.path import basename
    # now you can call it directly with basename
    #print basename(af)
    #str(af[:-4])+
    global theInputFiles

    from collections import OrderedDict
    
    print "choice 1 = ", af
    
    print "dict ipr2go = " , kf
    
    cfile = str(af[:-4])+"_counts_ipr.txt"
    cf = open(cfile, "w")
    
    iprs1 = []
    iprs2 = []
    gos1 = []
    gos2 = []
    iprn_dict = {}
    gon_dict = {}

    

 
    FileCheck(kf)
    kfile = open(kf,"r")
    fields = kfile.readline()
    for i in kfile.readlines():
        
        i = i.strip()
        entry = i.split("\t")
        if len(entry) > 1:
            nameipr = entry[1]
            namego = entry[-2]
            for j in entry:
                if j[:3] == "IPR":
                    #print j
                    if len(nameipr) > 50:
                        nameipr = nameipr[:50]
                    iprn_dict[j] = nameipr # [1]
                    
                elif j[:3] == "GO:":
                    #print j
                    if len(namego) > 50:
                        namego = namego[:50]
                    gon_dict[j] = namego
                    


    for i in kfile.readlines():
        i = i.strip()
        entry = i.split("\t")
        #print entry
        if len(entry) > 3 :
            namego = entry[-2]
            #print namego
            for j in entry:
                if j[:3] == "GO:":
                    #print j
                    gon_dict[j] = namego
    
    acs1 = []

    iacs1 = []
    iprc1 = 0
    goc1 = 0
    iprs11 = []
    gos11 = []
    all_ac1 = []
    all_ac2 = []
    gacs1 = []
    gacs2 = []
    notfound = 0
    num_lines1 = sum(1 for line in open(af))
    with open(af,"r") as f1:
        for i in f1.readlines()[1:]:
            #print i
            i = i.strip()
            items = i.split("\t")
            ac = items[0]
            if ac not in all_ac1:
                all_ac1.append(ac)
            for j in items:
                #print j
                if j[0:3] == "IPR":
                    iprc1 += 1
                    if ac not in iacs1:
                        iacs1.append(ac)                    
                    #print index(j)
                    #print j.split("; ")
                    for i in j.split("; "):
                        iprs1.append(i)
                        if i not in iprs11:
                            iprs11.append(i)
                elif j[0:3] == "GO:":
                    goc1 += 1
                    #print j.split("; ")
                    if ac not in gacs1:
                        gacs1.append(ac)
                    for i in j.split("; "):
                        gos1.append(i)
                        if i not in gos11:
                            gos11.append(i)
                #else:
                #notfound += 1
                
                

    iacs2 = []
    #acs2 = []
    iprc2 = 0
    goc2 = 0
    iprs22 = []
    gos22 = []

    no_go = len(all_ac1) - len(gacs1) 
    no_ipr = len(all_ac1) - len(iacs1) 
    

    #    PIE CHART


    
    
    
    c  = Counter(iprs1)
    print "\nTOTAL iprs1 found : ", len(iprs1) # 6597
    print "UNIQUE iprs1 : ", len(c), "\n", c.most_common(10) # 2088

    cf.write("########################## UNIQUE IPRS: "+str(len(c))+"\n")

    c_sorted_by_value = OrderedDict(sorted(c.items(), key=lambda x: x[1]), reverse = True)
    for i in c_sorted_by_value:
        if i[:3] == "IPR":
            if i in iprn_dict.keys():
                cf.write(str(i)+"\t"+str(c[i])+"\t"+str(iprn_dict[i])+"\n")
            else:
                cf.write(str(i)+"\t"+str(c[i])+"\n")
            
            # print i, f[i]
        
    #for i in dict(c.most_common()):
     #   cf.write(str(i)+"\t"+str(c[i])+"\n")
        # print i, c[i]
        
    ### COVERAGE IPR 1
    #print iprc1, acs1
    try:
        covi =  (len(iacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enri = (len(iprs11)/float(len(iprs1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation IPR 1 :\n\tAc Coverage: ", ("{0:.2%}".format(covi))#
        print "\tEnrichment", ("{0:.2%}".format(enri))
    except:
        print "No IPRs 1 found"



    # GO SUM 
    # COVERAGE GO 1

    
    f  = Counter(gos1)
    print "\nTOTAL GO IDS 1 found : ", len(gos1) # 6597
    print "UNIQUE GO IDs 1: ", len(f), "\n", f.most_common(10) # 1183

    cf.write("\n\n######################## UNIQUE GOS: "+str(len(f))+"\n")

    f_sorted_by_value = OrderedDict(sorted(f.items(), key=lambda x: x[1]), reverse = True)
    for i in f_sorted_by_value:
        if i[:3] == "GO:":
            if i in gon_dict.keys():
                cf.write(str(i)+"\t"+str(f[i])+"\t"+str(gon_dict[i])+"\n")
            else:
                cf.write(str(i)+"\t"+str(f[i])+"\n")
                # print i, f[i]

    try:
        covg = (len(gacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enrg = (len(gos11)/float(len(gos1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 1 :\n\tAc Coverage : ", ("{0:.2%}".format(covg))#
        print "\tEnrichment", ("{0:.2%}".format(enrg))
    except:
        print "No GOs 1 found!"

        
    #pie_plot_legends(ipr_pie_names, ipr_pie_nums)
    #pie_plot_legends(go_pie_names, go_pie_nums)
    # The slices will be ordered and plotted counter-clockwise.

    theInputFiles.append(cf.name)

    go_pie_names = []
    go_pie_names2 = []
    go_pie_nums = []
    
    fout = 0
    repercentgo = 0

    for i, j in f.most_common():
        if i[:3] == "GO:":
            go_pie_names.append(i)
            go_pie_nums.append(j)
            new = str(j) 
            if i in gon_dict.keys():
                if gon_dict[i] == True:
                    new = str(j) +" "+ gon_dict[i]
                else:
                    new = str(j) 
            go_pie_names2.append(new) # TO LABELS
            fout = fout + j 
            #print fout, len(iprs1), repercent
            repercentgo = fout / len(gos1)
            #half = int(len(gos1)*0.3)
            half = int(len(gos1)-fout)
            if len(go_pie_names) > 13 :
                go_pie_names.append("Other")
                go_pie_names.append("No GO ID")
                go_pie_nums.append(half)
                go_pie_nums.append(no_go)
                new = str(half) +" Other"
                new2 = str(no_go) + " No GO ID"
                go_pie_names2.append(new) # TO LABELS
                go_pie_names2.append(new2) # TO LABELS
                break


    repercent2go = fout / (len(gos11)+1)
    #print repercentgo, repercent2go
        

    

    ipr_pie_names = []
    ipr_pie_nums = []

    ipr_pie_names2 = []

    fout = 0
    repercentipr = 0
    for i, j in c.most_common():
        if i[:3] == "IPR":
            ipr_pie_names.append(i)
            ipr_pie_nums.append(j)
            if i in iprn_dict.keys():
                if iprn_dict[i] == True:
                    new = str(j) +" "+ iprn_dict[i]
                else:
                    new = str(j) 
            new = str(j) +" "+ iprn_dict[i]
            ipr_pie_names2.append(new) # TO LABELS
            fout = fout + j 
            #print fout, len(iprs1), repercent
            repercentipr = fout / len(iprs1)
            half = int(len(gos1)-fout)
            
            if len(ipr_pie_names) > 13 :
                #half = int(len(iprs1)-fout)
                ipr_pie_names.append("Other")
                ipr_pie_names.append("No IPR ID")
                
                ipr_pie_nums.append(half)
                
                
                new = str(half) +" Other"
                ipr_pie_names2.append(new) # TO LABELS

                if no_ipr<0:# no use
                    ipr_pie_nums.append(no_ipr)
                    new2 = str(no_ipr) + " No IPR ID"
                    ipr_pie_names2.append(new2) # TO LABELS
                break
            
            
    repercent2ipr = fout / (len(iprs11)+1)
    #print repercentipr, repercent2ipr

    fig = plt.figure()
    fig.canvas.set_window_title(str(cf.name[:-4]))

    
    ## IPR ##
    labels = ipr_pie_names # 'Frogs', 'Hogs', 'Dogs', 'Logs'
    labels2 = ipr_pie_names2
    sizes = ipr_pie_nums # [15, 30, 45, 10]

    ## GO ##
    #labels = go_pie_names
    #labels2 = go_pie_names2
    #sizes = go_pie_nums 
    
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    #explode = (0.1, 0, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90) #labels=labels
    # Set aspect ratio to be equal so that pie is drawn as a circle.

    plt.axis('equal')

    #plt.title(str(af))
    
    plt.legend(labels2, loc='best', shadow=True)
    pietitle = " ACsIPR= %s , TotIPRs= %s , UniqIPRs = %s"  % (len(iacs1),len(iprs1),len(iprs11))
    pietitle2 = "\nACs= %s, Coverage = %s , Enrichment = %s " % (len(all_ac1), "{0:.2%}".format(covi), "{0:.2%}".format(enri))

    
    plt.suptitle("20 Most Common IPR domains in Protein Accessions "+ pietitle+pietitle2)
    #plt.legend(title="técnica")

    #ax.legend([extra, bar_0_10, bar_10_100], ("My explanatory text", "0-10", "10-100"))

    #plt.savefig(name)
    
    repercentgo = repercentgo * 100
    centgo = "{0:.1f}".format(repercentgo)
    
    repercentipr = repercentipr * 100
    centipr = "{0:.1f}".format(repercentipr)

    
    
    #print "\nRepresenting %s %% of the total IPRs " % centipr
    #print "Representing %s %% of the total GOs " % centgo

    """from colour import Color
    red = Color("red")
    colors = list(red.range_to(Color("green"),10))
    """

    
    print"\nTOTAL UNIQUE ACS 1 :",len(all_ac1)
    
    print "\nAC with IPR : ", len(iacs1)
    print "TOTAL IPRs1 : ",len(iprs1)
    print "UNIQUE IPRs1 : ",len(iprs11)

    print "AC with GO : ", len(gacs1)
    print "TOTAL GOs1 : ",len(gos1)
    print "UNIQUE GOs1 : ",len(gos11)


    

    

    plt.show()

    
    

    print "(complete!)"

    """## GO ##

    labels2 = go_pie_names # 'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes2 = go_pie_nums # [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    #explode = (0, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes2, labels=labels2, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90)
    # Set aspect ratio to be equal so that pie is drawn as a circle.

    plt.axis('equal')

    plt.show()"""
        

    print "\n", cf.name, " saved!"

    
    
    #print ef, " saved!"
    #ef.close()
    #cff.close()
    #dff.close()
    kfile.close()
    cf.close()

#count_ipr_pie("output/tax5671_bmq_ac_ipr.txt", "output/_bmq_ipr_go.txt")
#250_down_all_filter_250_tax5671_bmq_ac_ipr_go
#count_ipr_pie("250_down_all_filter_250_tax5671_bmq_ac_ipr_go.txt", "output/_bmq_ipr_go.txt")

######################################################################################
kkf = "kegg_lin.txt"
def count_path_pie(af, kkf): # n = 1000

    """Using the most common dict (COLLECTIONS) file from last functions,
    we are able to build a default proteome set (bf) 
    and compare it with input proteome set (af) 
    combining de ipr2go EBI dictionary (kf)                                     """


    #from os.path import basename
    # now you can call it directly with basename
    #print basename(af)
    #str(af[:-4])+
    
    global theInputFiles

    from collections import OrderedDict

    numbers = ["1","2","3","4","5","6","7","8","9","0"]
    minletters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","k","r","s","t","u","v","x","y","z"]
    
    
    print "choice 1 = ", af
    
    print "dict ipt2go = " , kf
    print "dict path = " , kkf
    
    cfile = str(af[:-4])+"_counts_paths_.txt"
    cf = open(cfile, "w")
    
    iprs1 = []
    iprs2 = []
    gos1 = []
    gos2 = []
    iprn_dict = {}
    gon_dict = {}
   

    FileCheck(kf)
    kfile = open(kf,"r")
    for i in kfile.readlines():
        i = i.strip()
        entry = i.split("\t")
        if len(entry) > 1:
            nameipr = entry[1]
            namego = entry[-2]
            for j in entry:
                if j[:3] == "IPR":
                    #print j
                    iprn_dict[j] = nameipr # [1]
                elif j[:3] == "GO:":
                    #print j
                    if len(namego) > 50:
                        namego = namego[:50]
                    gon_dict[j] = namego
                #elif j[8] == ":" and j[7] in numbers and j[1] in minletters:
                    
                    


    for i in kfile.readlines():
        i = i.strip()
        entry = i.split("\t")
        #print entry
        if len(entry) > 3 :
            namego = entry[-2]
            #print namego
            for j in entry:
                if j[:3] == "GO:":
                    #print j
                    gon_dict[j] = namego
    
    acs1 = []

    iacs1 = []
    iprc1 = 0
    goc1 = 0
    iprs11 = []
    gos11 = []
    all_ac1 = []
    all_ac2 = []
    gacs1 = []
    gacs2 = []
    notfound = 0
    num_lines1 = sum(1 for line in open(af))

    with open(af,"r") as f1:
        for i in f1.readlines()[1:]:
            #print i
            i = i.strip()
            items = i.split("\t")
            ac = items[0]
            if ac not in all_ac1:
                all_ac1.append(ac)
            for j in items:
                #print j
                if j[0:3] == "IPR":
                    iprc1 += 1
                    if ac not in iacs1:
                        iacs1.append(ac)                    
                    #print index(j)
                    #print j.split("; ")
                    for i in j.split("; "):
                        iprs1.append(i)
                        if i not in iprs11:
                            iprs11.append(i)
                elif j[0:3] == "GO:":
                    goc1 += 1
                    #print j.split("; ")
                    if ac not in gacs1:
                        gacs1.append(ac)
                    for i in j.split("; "):
                        gos1.append(i)
                        if i not in gos11:
                            gos11.append(i)
                #else:
                #notfound += 1
                
                

    iacs2 = []
    #acs2 = []
    iprc2 = 0
    goc2 = 0
    iprs22 = []
    gos22 = []

    no_go = len(all_ac1) - len(gacs1)
    no_ipr = len(all_ac1) - len(iacs1)

    #print no_go, no_ipr
    #print len(gacs1), len(iacs1)
    #print 
    
                

    

    ##############3   PIE CHART     #######################3
    
    
    
    c  = Counter(iprs1)
    print "\nTOTAL iprs1 found : ", len(iprs1) # 6597
    print "UNIQUE iprs1 : ", len(c), "\n", c.most_common(10) # 2088

    cf.write("########################## UNIQUE IPRS: "+str(len(c))+"\n")

    c_sorted_by_value = OrderedDict(sorted(c.items(), key=lambda x: x[1]), reverse = True)
    for i in c_sorted_by_value:
        if i[:3] == "IPR":
            if i in iprn_dict.keys():
                cf.write(str(i)+"\t"+str(c[i])+"\t"+str(iprn_dict[i])+"\n")
            else:
                cf.write(str(i)+"\t"+str(c[i])+"\n")            

            
            # print i, f[i]
        
    #for i in dict(c.most_common()):
     #   cf.write(str(i)+"\t"+str(c[i])+"\n")
        # print i, c[i]
        
    ### COVERAGE IPR 1
    #print iprc1, acs1
    try:
        cov =  (len(iacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enr = (len(iprs11)/float(len(iprs1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation IPR 1 :\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No IPRs 1 found"



    ### GO SUM ###
    ### COVERAGE GO 1

   
    
    f  = Counter(gos1)
    print "\nTOTAL GO IDS 1 found : ", len(gos1) # 6597
    print "UNIQUE GO IDs 1: ", len(f), "\n", f.most_common(10) # 1183

    cf.write("\n\n######################## UNIQUE GOS: "+str(len(f))+"\n")

    f_sorted_by_value = OrderedDict(sorted(f.items(), key=lambda x: x[1]), reverse=True)
    for i in f_sorted_by_value:
        if i[:3] == "GO:":
            if i in gon_dict.keys():
                cf.write(str(i)+"\t"+str(f[i])+"\t"+str(gon_dict[i])+"\n")
            else:
                cf.write(str(i)+"\t"+str(f[i])+"\n")
                # print i, f[i]

    try:
        cov = (len(gacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enr = (len(gos11)/float(len(gos1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 1 :\n\tAc Coverage : ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No GOs 1 found!"

        
    #pie_plot_legends(ipr_pie_names, ipr_pie_nums)
    #pie_plot_legends(go_pie_names, go_pie_nums)
    # The slices will be ordered and plotted counter-clockwise.

    theInputFiles.append(cf.name)

    go_pie_names = []
    go_pie_names2 = []
    go_pie_nums = []
    
    fout = 0
    repercentgo = 0

    for i, j in f.most_common():
        if i[:3] == "GO:":
            go_pie_names.append(i)
            go_pie_nums.append(j)
            new = str(j) 
            if i in gon_dict.keys():
                new = str(j) +" "+ gon_dict[i]
            
            go_pie_names2.append(new) # TO LABELS
            fout = fout + j
            #print fout, len(iprs1), repercent
            repercentgo = fout / len(gos1)
            
            #half = int(len(gos1)*0.3)
            half = int(len(gos1)-fout)
            if len(go_pie_names) > 20 :
                go_pie_names.append("Other")
                go_pie_names.append("No GO ID")
                go_pie_nums.append(half)
                go_pie_nums.append(no_go)
                new = str(half) +" Other"
                new2 = str(no_go) + " No GO ID"
                go_pie_names2.append(new) # TO LABELS
                go_pie_names2.append(new2) # TO LABELS
                break


    repercent2go = fout / (len(gos11)+1)
    print repercentgo, repercent2go
        

    

    ipr_pie_names = []
    ipr_pie_nums = []

    ipr_pie_names2 = []

    fout = 0
    repercentipr = 0
    for i, j in c.most_common():
        if i[:3] == "IPR":
            ipr_pie_names.append(i)
            ipr_pie_nums.append(j)
            new = str(j) +" "+ iprn_dict[i]
            ipr_pie_names2.append(new) # TO LABELS
            fout = fout + j
            #print fout, len(iprs1), repercent
            repercentipr = fout / len(iprs11)
            if repercentipr > 0.5 :
                break
            
    repercent2ipr = fout / (len(iprs11)+1)
    #print repercentipr, repercent2ipr

    fig = plt.figure()
    fig.canvas.set_window_title(str(cf.name[:-4])) # no .name

    
    ## IPR ##
    #labels = ipr_pie_names # 'Frogs', 'Hogs', 'Dogs', 'Logs'
    #labels2 = ipr_pie_names2
    #sizes = ipr_pie_nums # [15, 30, 45, 10]

    ## GO ##
    labels = go_pie_names
    labels2 = go_pie_names2
    sizes = go_pie_nums 
    
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    explode = (2, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90)
    # Set aspect ratio to be equal so that pie is drawn as a circle.

    plt.axis('equal')

    #plt.title(str(af))
    
    plt.legend(labels2, loc='best', shadow=True)
    #pietitle = "ACsGO=" +str(len(gacs1))+", TotGOs="+str(len(gos1))
    pietitle = " ACsGO= %s , TotGOs= %s , UniqGOs = %s"  % (len(gacs1),len(gos1), len(gos11))
    
    plt.suptitle("20 Most Common GO IDs in Protein Accessions " + pietitle)
    #plt.legend(title="técnica")

    #ax.legend([extra, bar_0_10, bar_10_100], ("My explanatory text", "0-10", "10-100"))

    #plt.savefig(name)
    #print "Representing %s of the total IPRs " % repercentipr

    repercentgo = repercentgo * 100
    centgo = "{0:.1f}".format(repercentgo)

    repercent2ipr= repercent2ipr * 100
    centipr = "{0:.1f}".format(repercent2ipr)

    
    
    #print "Representing %s %% of the total IPRs " % centipr
    #print "Representing %s %% of the total GOs " % centgo

   

    plt.show()

    
    

    print "(complete!)"

    """## GO ##

    labels2 = go_pie_names # 'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes2 = go_pie_nums # [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    #explode = (0, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes2, labels=labels2, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90)
    # Set aspect ratio to be equal so that pie is drawn as a circle.

    plt.axis('equal')

    plt.show()"""

    print "\nAC with IPR : ", len(iacs1)
    print "TOTAL IPRs1 : ",len(iprs1)
    print "UNIQUE IPRs1 : ",len(iprs11)
    

    print "AC with GO : ", len(gacs1)
    print "TOTAL GOs1 : ",len(gos1)
    print "UNIQUE GOs1 : ",len(gos11)


    print"\nTOTAL UNIQUE ACS 1 :",len(all_ac1)
    print"\nTOTAL UNIQUE ACS 1 :",len(all_ac2)

        

    print "\n", cf.name, " saved!"

    
    
    #print ef, " saved!"
    #ef.close()
    #cff.close()
    #dff.close()
    kfile.close()
    cf.close()


#count_go_pie("250_down_all_filter_250_tax5671_bmq_ac_ipr_go.txt", "output/_bmq_ipr_go.txt")


######################################################################################
#### BAR PLOT
######################################################################################

def from_file_get_uniques_bar(af):

    # Gets tax id and retrieves uni ids list from uniprot web
    # then enrichment analysis with list of accessions

    #from os.path import basename
    # now you can call it directly with basename
    #print basename(af)
    #str(af[:-4])+
    global theOutputbars
    

    from collections import OrderedDict
    
    print "\n########################################\nChoice 1 = ", af
    
    #print "dict ipr2go = " , kf
    
    mfile = str(af[:-4])+"_ipr_go_counts.txt"
    mf = open(mfile, "w+")
    
    iprs1 = []
    iprs2 = []
    gos1 = []
    gos2 = []
    iprn_dict = {}
    gon_dict = {}

    
    
    acs1 = []

    iacs1 = []
    iprc1 = 0
    goc1 = 0
    iprs11 = []
    gos11 = []
    all_ac1 = []
    all_ac2 = []
    gacs1 = []
    gacs2 = []
    # num_lines1 = sum(1 for line in open(af))
    with open(af,"r") as f1:
        for i in f1.readlines():
            #print i
            i = i.strip()
            items = i.split("\t")
            ac = items[0]
            #acs1.append(ac)
            if ac not in all_ac1:
                all_ac1.append(ac)
            for j in items:
                #print j
                
                if j[0:3] == "IPR":
                    iprc1 += 1
                    if ac not in iacs1:
                        iacs1.append(ac)                    
                    #print index(j)
                    #print j.split("; ")
                    for i in j.split("; "):
                        iprs1.append(i)
                        if i not in iprs11:
                            iprs11.append(i)
                elif j[0:3] == "GO:":
                    goc1 += 1
                    #print j.split("; ")
                    if ac not in gacs1:
                        gacs1.append(ac)
                    for i in j.split("; "):
                        gos1.append(i)
                        if i not in gos11:
                            gos11.append(i)
                
                

    
    
    
                
    ##############3   BAR CHART     #######################3

    print "TOTAL IPRs : ",len(iprs1)
    print "TOTAL GOs : ",len(gos1)

    print "UNIQUE IPRs : ",len(iprs11)
    print "UNIQUE GOs : ",len(gos11)


    print"\nTOTAL ACS  :",len(all_ac1)
    print "\nAC  with IPR :",len(iacs1)     # Unique ACS, that have ipr
    print "\nAC  WITH GO :",len(gacs1) 

    
    
    
    c  = Counter(iprs1)
    print "\nTOTAL IPRS  found : ", len(iprs1) # 6597
    print "UNIQUE IRPS : ", len(c), "\n", c.most_common(10) # 2088

    mf.write("########################## UNIQUE IPRS: "+str(len(c))+"\n")

    c_sorted_by_value = OrderedDict(sorted(c.items(), key=lambda x: x[1]))
    for i in c_sorted_by_value:
        if i[:3] == "IPR":
            if i in iprn_dict.keys():
                mf.write(str(i)+"\t"+str(c[i])+"\t"+str(iprn_dict[i])+"\n")
                # print i, f[i]
            else:
                mf.write(str(i)+"\t"+str(c[i])+"\n")
                # print i, f[i]
        
    #for i in dict(c.most_common()):
     #   cf.write(str(i)+"\t"+str(c[i])+"\n")
        # print i, c[i]
        
    ### COVERAGE IPR 1
    #print iprc1, acs1
    try:
        cov =  (len(iacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enr = (len(iprs11)/float(len(iprs1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation IPR 1 :\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No IPRs 1 found"



    ### GO SUM ###
    ### COVERAGE GO 1

   
    
    f  = Counter(gos1)
    print "\nTOTAL GO IDS 1 found : ", len(gos1) # 6597
    print "UNIQUE GO IDs 1: ", len(f), "\n", f.most_common(10) # 1183

    mf.write("\n\n######################## UNIQUE GOS: "+str(len(f))+"\n")

    f_sorted_by_value = OrderedDict(sorted(f.items(), key=lambda x: x[1]))
    for i in f_sorted_by_value:
        if i[:3] == "GO:":
            if i in gon_dict.keys():
                mf.write(str(i)+"\t"+str(f[i])+"\t"+str(gon_dict[i])+"\n")
            else:
                mf.write(str(i)+"\t"+str(f[i])+"\n")
                # print i, f[i]

    try:
        cov = (len(gacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enr = (len(gos11)/float(len(gos1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 1 :\n\tAc Coverage : ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No GOs 1 found!"

        
    #pie_plot_legends(ipr_pie_names, ipr_pie_nums)
    #pie_plot_legends(go_pie_names, go_pie_nums)
    # The slices will be ordered and plotted counter-clockwise
    
    
    return all_ac1, iacs1, gacs1, iprs1, iprs11, gos1, gos11, mf.name

    

###########################################################################3

def bar_plot2labels(all_files):


    #!/usr/bin/env python
    #http://matplotlib.org/examples/api/barchart_demo.html
    global theInputFiles

    col_ac = ()
    col_iac = ()
    col_gac = ()

    col_ipr = ()
    col_uniq_ipr = ()
    col_go = ()
    col_uniq_go = ()
    
  

    bar_names = ()
    out_files = []
     

    for i in all_files:
        all_ac, iacs, gacs, iprs, iprs_q, gos, gos_q, outf = from_file_get_uniques_bar(i)
        
        col_ac = col_ac + (len(all_ac),)
        col_iac = col_iac + (len(iacs),)
        col_gac = col_gac + (len(gacs),)
        
        col_ipr = col_ipr + (len(iprs),)
        col_uniq_ipr = col_uniq_ipr + (len(iprs_q),)
        
        col_go = col_go + (len(gos),)
        col_uniq_go = col_uniq_go + (len(gos_q),)
        out_files.append(outf)

        # Create Title Name

        if "/" in i:
            i = i.split("/")
            name1 = i[1]
            bar_name = name1.split("_")[0]
            bar_names = bar_names + (bar_name,)
        elif "\\" in i:
            i = i.split("\\")
            name1 = i[1]
            bar_name = name1.split("_")[0]
            bar_names = bar_names + (bar_name,)
        else:
            i = i.split("\\")
            name1 = i[0]
            bar_name = name1.split("_")[0]
            bar_names = bar_names + (bar_name,)
                        
            
            
        


    # a bar plot with errorbars
    N = len(all_files)

    print col_uniq_ipr, col_uniq_go


    menStd =   (2, 3, 4, 1, 2)

    ind = np.arange(N)  # the x locations for the groups

    width = 0.2       # the width of the bars

    #fig = plt.figure()
    


    fig, ax = plt.subplots()
    fig.canvas.set_window_title("Bar Plot")
    
    #rects1 = ax.bar(ind, menMeans, width, color='r', yerr=menStd)


    rects01 = ax.bar(ind-width, col_ac, width, color='w' ) #ind, , width, color='w')
    rects02 = ax.bar(ind-width,  col_iac, width, color='m') # ind,, width, color='m')
    rects03 = ax.bar(ind-width, col_gac, width, color='c') # ,align='center') # ind, col_gac, width, color='c')
    

    rects1 = ax.bar(ind, col_ipr, width, color='g')
    rects3 = ax.bar(ind, col_uniq_ipr, width, color='r')
 

    womenStd =   (3, 5, 2, 3, 3)
    #rects2 = ax.bar(ind+width, womenMeans, width, color='y', yerr=womenStd)

    
    rects2 = ax.bar(ind+width, col_go, width, color='y')
    rects4 = ax.bar(ind+width, col_uniq_go, width, color='b')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Number of IDs')
    ax.set_title('Scores by file and type')
    ax.set_xticks(ind)
    ax.set_xticklabels( bar_names) #  rotation=45) # most 4 args


    print rects01[0], rects02[0], rects03[0], rects2[0]

    ax.legend( (rects01[0], rects02[0], rects03[0],rects1[0], rects3[0], rects2[0], rects4[0]) ,  ('Uniq_AC', 'AC IPR', 'AC GO','IPR', 'Uniq_IPR', 'GO','Uniq_GO') )

    #autolabel(rects1)
    # attach some text labels
        
    for rect in rects01:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')

    for rect in rects02:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')
        
    for rect in rects03:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')


    
    for rect in rects1:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')

    for rect in rects3:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')


    #autolabel(rects2)

    # attach some text labels

    for rect in rects2:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')
        
    for rect in rects4:
        height = rect.get_height()
        ax.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')

    onetitle = "\t".join(all_files)

    plt.title(str(onetitle))
    for i in out_files:
        print i, " saved!"
        theInputFiles.append(i)

    print "(complete!)"
    
    plt.show()
    # plt.savefig("figure.pdf")

    

    """fig = plt.figure()
    width = .35
    ind = np.arange(len(OY))
    plt.bar(ind, OY)
    plt.xticks(ind + width / 2, OX)

    fig.autofmt_xdate()

    plt.savefig("figure.pdf")"""
    
    


#af,bf,cf,df,ef, kf = "output/tax5671_bs_uniprot1.txt","output/tax5671_bmq_ac_ipr_go.txt", "output/tax5664_bs_uniprot.txt","output/tax5661_bs_uniprot.txt","output/tax5664_bs_uniprot.txt","output/_bmq_ipr_go.txt"

#all_files = "output/tax5671_bs_uniprot.txt","output/tax5671_bmq_ac_ipr_go.txt", "output/tax5664_bs_uniprot.txt","output/tax5661_bs_uniprot.txt"# ,"output/tax5664_bs_uniprot.txt","output/_bmq_ipr_go.txt"

#files = ["output/tax5671_bs_uniprot1.txt","output/tax5661_bs_uniprot1.txt","output/tax5664_bs_uniprot1.txt","output/tax5671_bmq_ac_ipr_go.txt","output/_bmq_ipr_go.txt"]


#############################
from bioservices import *

def get_bmq_ipr2go():

    global theInputFiles

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
    if bFile.name not in theInputFiles:
        theInputFiles.append(bFile.name)

    print bFile, " saved!"
    
    
def FileCheck(kf):
    try:
      open(kf, "r")
      return 1
    except IOError:
      print "Error: File %s not appear to exist.\nLet 's create <_bmq_ipr_go.txt>" % kf
      get_bmq_ipr2go()
      return 0

    
################33##########3

#from bioservices_functions_tk import FileCheck
def from_go_get_names(go_ids):

   import collections

   iprs1 = []
   iprs2 = []
   gos1 = []
   gos2 = []
   iprn_dict = {}
   gon_dict = {}
   go_names = {}

   kf = "output/_bmq_ipr_go.txt"

   FileCheck(kf)
   kfile = open(kf,"r")
   for i in kfile.readlines():
      i = i.strip()
      entry = i.split("\t")
      if len(entry) > 1:
         nameipr = entry[1]
         namego = entry[-2]
         for j in entry:
            if j[:3] == "IPR":
               #print j
               iprn_dict[j] = nameipr # [1]
            elif j[:3] == "GO:":
               #print j
               gon_dict[j] = namego

   for i in go_ids:
      if i in gon_dict.keys():
         go_names[i] = gon_dict[i]
      else:
         go_names[i] = "null"
         
         
      
   
   
   return go_names



###############################################################################
#####   SCATTER
###############################################################################
def scatter_ipr_functions(sf):

   print "\nWork in process . . . "


def scatter_path(sf):

   print "\nWork in process . . . "
   

def scatter_go_functions(sf):

   # http://thomas-cokelaer.info/blog/2012/04/481/

    from collections import OrderedDict
    from collections import defaultdict
    from collections import Counter
    import numpy as np
    
    

    from operator import itemgetter

    global theInputFiles
    #time_begin = datetime.datetime.fromtimestamp(time.time())

    

    # Open data to plot

    dfile = open(sf, "r")
    first = dfile.readline()
    first = first.strip()
    fields = first.split("\t")
    print fields
    positions = ["p-value","odds_ratio","benjamini", "bonferroni"]
    for i in fields:
        if "p-value" == i:
            ind_pv = fields.index("p-value")
            print i
        elif "odds_ratio" == i:
            ind_or = fields.index("odds_ratio")
            print i
        elif "Medium/Light" == i:
            ind_ml = fields.index("Medium/Light")
            print i
        elif "Area" == i:
            ind_ar = fields.index("Area")
            print i
        else:
            pass
            

    pvals6 = []
    counts6 = []
    bjmn6 = []
    
    bfrn6 = []

    
    
    gos = []
    counts = []
    pvals = []
    lines = []
    lines_id = {}
    count2pv_dict = defaultdict(list)
    dict_ipr2pv = defaultdict(list)
    dict_ipr2cts = defaultdict(list)
    dict_ipr2or = defaultdict(list)
    areas = []
    medlights =[]
    ac_ids = []
    dict_ac2ml= {}
    dict_ac2ar= {}
    dict_ipr2area = {}
    dict_ac2lig = {}
    go_ids = []
    gogo_ids = []
    dict_ac2go = defaultdict(list)
    dict_ac2line = {}
    for i in dfile.readlines():
        i = i.strip()
        #print i
        #if "p-value" in fields: # IPR NAME
        #    pos = fields.index('p-value')
        #print pos
        #    ipr_name = line[pos] if len(line) > 2 else "null"
        #    ipr_names1.append(ipr_name)
        items = i.split("\t")
        
        
        #lines.append(items)
        # l = [[0, 1, 'f'], [4, 2, 't'], [9, 4, 'afsd']]
        # l.sort(key=lambda x: x[2])
        if len(items) > 8:
            #print items
            ac = items[0]
            ac_ids.append(ac)
            #gos.append(items[0])
            lines.append(items)
            ipr = items[0]
            #pv =items[ind_pv]
            #odds = i[ind_or]
            medlig = items[ind_ml]
            area = items[ind_ar]
            #print area, medlig
            
            valy = str(area)
            y2  = valy.replace(",",".")
            
            if y2:
               y = np.log2(float(y2))
               areas.append(y)
            else:
               y = 0
               areas.append(y)
            
            valx = str(medlig)
            x2  = valx.replace(",",".")
            #print x2

            
            if x2:
               x = np.log2(float(x2)) 
               medlights.append(x)
            else:
               x = 0
               medlights.append(x)
               
            
            #pvals.append(items[-2])
            cts = items[2]
            #print cts, pv
            #pvals6.append(pv)
            #counts6.append(cts)
            dict_ac2lig[ac] = x
            
           
            #print cts
            #counts.append(items[2])
            #count2pv_dict[cts].append(pv) #= pv
            #count2pv_dict[pv] = cts# append(pv)   USING DICT of P_VALUES
            joint = "\t".join(items)
            
            lines_id[ipr]= joint
            #print joint
            #dict_ipr2pv[ipr] = str(pv)
            #dict_ipr2cts[ipr] = str(cts)
            #dict_ipr2or[ipr] = str(odds)

            dict_ac2ml[ac] = x
            dict_ac2ar[ac] = y
            dict_ac2line[ac] = joint

            for i in items:
               if i[:3] == "GO:":
                  gos = i.split("; ")
                  for i in gos:
                      dict_ac2go[ac].append(i)
                      #if i not in go_ids:
                      go_ids.append(i)
                      if i not in gogo_ids:
                          gogo_ids.append(i)
                  
                  
           
      
            
        else:
            print len(items), items
        
            

    #floats1 = " ".join(pvals)
    #floats = map(float, floats1.split())

    #print len(counts)
    #print len(pvals)
    
    #pvals = sorted(pvals)
    
    #counts = sorted(counts)
    
    bjmn = []
    
    bfrn = []
    d_sorted_by_value = OrderedDict(sorted(count2pv_dict.items(), key=lambda x: x[0]))
    for k, v in d_sorted_by_value.items():
        #print "%s: %s" % (k, v)
        pvals.append(k)
        counts.append(v)
        
           

    #pv_sorted = OrderedDict(sorted(dict_ipr2pv.itens(), key=lambda: x: x[1]))
    #it = iter(sorted(d.iteritems()))
    #for key in sorted(dict_ipr2pv):
    #    print dict_ipr2pv[key]

    
    
    pvals2 = []
    counts2 = []
    
    lf = open(str(dfile.name[:-4])+"_order_scatter.txt","w")
    header1 = "\t".join(fields)
    lf.write(header1+"\n")

    #mf = open(str(dfile.name[:-4])+"_order_or.txt","w")
    #mf.write(header1+"\n")

    # import operator
    # x = {1: 2, 3: 4, 4: 3, 2: 1, 0: 0}
    # sorted_x = sorted(x.items(), key=operator.itemgetter(0))
    # sorted_x = sorted(x.items(), key=operator.itemgetter(1))
    # for w in sorted(d, key=d.get, reverse=True):
    # print w, d[w]
    # sorted(d.items(), key=lambda x: x[1])


    pvals5 = []
    counts5 = []


  
    

    
    order_by_ml = sorted(lines, key=itemgetter(ind_ml)) #lines.sort(key=lambda x: x[-2]) :error
    print len(order_by_ml), len(lines)
    
    pvals3 = []
    counts3 = []
    #pvs_sorted = OrderedDict(sorted(dict_ipr2pv.items(), key=lambda x: x[1]))

 
        

        
    # CORRECT VERSION
    bjmn2 = []
    bfrn2 = []
    medlig2 = []
    areas2 = []
    for i in order_by_ml: # list inside list
        #print i
        if len(i) > 8:
            val = i[ind_ml]  #else "999"
            #print val
            new_line2 = "\t".join(i)
            lf.write(new_line2+"\n")
            #pvals2.append(i[ind_pv])
            #counts2.append(i[2])
            medlig2.append(i[ind_ml])
            areas2.append(i[ind_ar])
            
                          
        else:
            print i

    lf.close()
    #mf.close()
    #print dfile.name[:-4]
    dfile.close()
    lines = []
    #global theInputFiles
    
    theInputFiles.append(lf.name) # PV - IS NOT IMPORTING FILE TO LIST
    #theInputFiles.append(mf.name) # OR
    #print theInputFiles
    #updatecombolist(theFilenames1)

    #dict_ac2ipr = from_ac_get_ipr(ac_ids)
    #dict_ipr2fam = from_ipr_get_tree(dict_ac2ipr)

    import matplotlib.pyplot as plt
    from numpy.random import random
    import random
    import matplotlib.cm as cm
    #plt.scatter(x, y, c=t, cmap=cm.colormap_name)

    
    #for i, j in dict_ipr2ml.iteritems():
    #lo = plt.scatter(i, j, marker='x', color=colors[0])
    
    one = Counter(medlights)
    medlights1 =[]
    areas1 = []
    acdata = []
    labels = []
    #colors1 = iter(cm.rainbow(np.linspace(0, 1, len(ys))))
    #for y in ys:
    #plt.scatter(x, y, color=next(colors))
    
    two = Counter(go_ids)
    count = 0
    print "TOTAL Unique IDs : ", len(two)
    colors = ['b', 'c', 'y', 'm', 'r','w','g']
    marks = ["x","o"]

    t = np.arange(10)
    #colors3 = cm.rainbow(np.linspace(0, 1, len(two)))
    #for y, c in zip(two, colors3):
       #plt.scatter(x, y, color=c)

    #colors = plt.cm.coolwarm(scaled_z)
    #import matplotlib.pyplot as plt
    cm = plt.cm.get_cmap('RdYlBu')
    #xy = range(20)
    #z = xy
    #sc = plt.scatter(xy, xy, c=z, vmin=0, vmax=20, s=35, cmap=cm)
    #plt.colorbar(sc)
    #plt.show()
    import matplotlib as mpl
    tot = 10
    
    #x, y, colors = np.random.random((3,10))
    #fig, ax = plt.subplots()
    #ax.scatter(x, y, c=colors, s=50, cmap=mpl.cm.Reds)
    mymap = plt.get_cmap("Reds")
    spaced_colors = linspace(0, 1, tot)
    #print spaced_colors
    import matplotlib.pyplot as mpl_plt
    import brewer2mpl

    # Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
    #set2 = brewer2mpl.get_map('Set2', 'qualitative', 8).mpl_colors
    set2 = brewer2mpl.get_map('Reds', 'sequential', 9).mpl_colors

    new_file = open(str(sf[:-4])+"_scatter.txt","w")

    fig = plt.figure()
    fig.canvas.set_window_title(str(new_file.name[:-4])) 
    

    #root3=Tk()
    #root3.title(str(dfile.name))

    
    colours = plt.cm.rainbow(np.linspace(0, 1, tot))
    #print colours
    mycolors = plt.cm.rainbow(np.linspace(0, 1, tot))
    
    

    dict_go2name = from_go_get_names(go_ids)
    
    #print two.most_common(1)
    power = []
    prev_j = 1
    for i, j in two.most_common(tot):
       print i, dict_go2name[i], j
       for k,v in dict_ac2go.iteritems():
          #print k, v
          for vv in v:
             if i == vv:
                #print k, vv, i, j, dict_ac2ml[k]
                if k not in acdata:
                    medlights1.append(dict_ac2ml[k])
                    areas1.append(dict_ac2ar[k])
                    acdata.append(k)
                    power.append(j)
                    new_file.write(dict_ac2line[k])

                #labels.append(i)
                #print labels
       #data = plt.cm.jet(data[medlights1, areas1])
       minn = 1
       maxn = 7
       nn = j /7
       #print nn
       nn = minn if nn < minn else maxn if nn > maxn else nn
       
       
       n  = j / nn
       #prev_j = j
       n = minn if n < minn else maxn if n > maxn else n
       #print set2[n], j, n, len(set2)
       print  acdata
       #top = max(medlights1)

       #colorme = [str(item/top) for item in medlights1]
       #plt.scatter(x, y, s=500, c=color) # http://stackoverflow.com/questions/8202605/matplotlib-scatterplot-colour-as-a-function-of-a-third-variable
       cc= mycolors[0]
       
       lk  = plt.scatter(medlights1, areas1, c = cc, marker = "o", s=50, label = str(dict_go2name[i]))# dict_go2name[i], color=random.choice(colors),

       
       #lk  = plt.scatter(medlights1, areas1, color = set2[n], marker = "o", s=50, label = str(i))# dict_go2name[i], color=random.choice(colors), 
       #plt.legend(lk, label=str(i),scatterpoints=1,loc='lower left', ncol=3,fontsize=8)
       #sc = plt.scatter(xy, xy, c=z, vmin=0, vmax=20, s=35, cmap=cm)
       #plt.colorbar(lk)
       
       #print len(medlights1),len(areas1)#,len(labels)
       medlights1 = []
       areas1 = []
       labels = []
       power = []
       acdata = []
       
       mycolors = np.delete(mycolors,0, 0 ) # np.delete(arr, [0,2,4], axis=0), 
       
    #
    
       
    #for i, j, k in zip(medlights1,areas1,labels):
       #lk  = plt.scatter(i, j, marker='o', color=colors[2], label='Low Outlier')
    
       

    #import operator
    #x = {1: 2, 3: 4, 4: 3, 2: 1, 0: 0}
    #sorted_x = sorted(one.items(), key=operator.itemgetter(0),reverse=True)
       
    #ll = plt.scatter(random(10), random(10), marker='o', color=colors[0], label='Low Outlier')
    #l  = plt.scatter(medlights1, areas1, marker='o', color=colors[1])
    #a  = plt.scatter(random(10), random(10), marker='o', color=colors[2])
    #h  = plt.scatter(random(10), random(10), marker='o', color=colors[3])
    #hh = plt.scatter(random(10), random(10), marker='o', color=colors[4])
    #ho = plt.scatter(random(10), random(10), marker='x', color=colors[4])

    #plt.legend((ll, a, h, hh, ho,lk),
       #        ('LoLo', 'Average', 'Hi', 'HiHi', 'High Outlier',"new"),
       #        scatterpoints=1,
       #        loc='lower left',
       #        ncol=3,
       #      fontsize=8)
         
    
    #plt.gray()
    #title = fields[0] +" vs. " +fields[ind_ml] + "\n" + dfile.name
    # create the general figure
    #fig1 = figure()
    #fig1.suptitle(str(title))

    
    #plt.figure(str(new_file.name))
    #fig = pylab.gcf()
    #fig.canvas.set_window_title('Test')
    
    plt.title('Scatter GO IDs')
    plt.ylabel('Log 2 (Area)')
    plt.xlabel('Log 2 (Medium/Light)')
    plt.legend(loc='upper left', scatterpoints=1, ncol=4, fontsize=8)

    # DRAW VERTIVAL LINE (70,100) to (70, 250)
    #ax1.plot([0, len(counts2)], [0.05, 0.05], 'g-', lw=2)
    
    ax1 = fig.add_subplot(111)  
    ax1.axvline(1, color='k', linestyle='--')

    ax2 = fig.add_subplot(111)  
    ax2.axvline(-1, color='k', linestyle='--')

    ax0 = fig.add_subplot(111)  
    ax0.axvline(0, color='k', linestyle='-')

    #fig.savefig(fig.my_figure_title + ".png")
    
    # draw vertical line from (70,100) to (70, 250)
    #ax1.plot([0, len(counts2)], [0.05, 0.05], 'g-', lw=2)
    
    plt.show()
    #root3.mainloop()
    

    print new_file.name," saved! "
    theInputFiles.append(new_file.name)

#sf16 = "250_up_all.txt"
#sf16 = "output/2peps_alls.txt"
#sf16 = "start\\2peps_log2dot.txt"



#scatter_go_functions(sf16)
# http://matplotlib.org/examples/color/colormaps_reference.html



   
time_end = datetime.datetime.fromtimestamp(time.time())
print("Time elapsed: ", str(time_end - time_begin))
