
import sys, os
from glob import *


import random
import winsound


import collections
from collections import Counter

import numbers
import decimal
import matplotlib as plt

from scipy import stats

#import tk_tax_id_v154

from bioservices_functions_tk import *

global theInputFiles

########################################################################

def merge_ac_with_ipr_ipr(file1,file2):

    global theInputFiles
    from collections import defaultdict

    info = "Select bigger file ( ComboBox 1) for yield matching. File 1 in Column 1 ."

    # METHOD SPECIFIC FOR BIOMART FORMAT, ONE IPR PER LINE

    uni_ids = []
    
    af = open(file1, "r") # choice 1
    bf = open(file2, "r") # choice 2

    ac1_dict = defaultdict(list)
    spacer = 77*"#"
    print spacer, "\n", info, "\n\n\tMerging ACIPR IDs with IPR IDS \n\t" ,file1 ," VS. ", file2
    
    acs1 = []
    fields =  af.readline()
    fields = fields.strip()

    iprs1_dict = defaultdict(list)
    iprs2_dict = defaultdict(list)
    iprs1 = []
    iprs2= []
    
    for i in af.readlines():
        i = i.strip()
        items = i.split("\t")
        ac1 = items[0]
        
        
        left = items[1:] if len(items)>0 else "null"
        ac1_dict[ac1].append(i)
        if ac1 not in acs1:
            acs1.append(ac1)
        for j in items:
            if j[:3] == "IPR":
                if i not in iprs1:
                    ipr = j
                    iprs1_dict[ipr].append(i)
                    iprs1.append(j)


    ac2_dict = defaultdict(list)
    acs2 = []
    fields2 =  bf.readline()
    fields2 = fields2.strip()
    
    print fields, "\t", fields2
    
    for i in bf.readlines()[1:]:
        i = i.strip()
        items = i.split("\t")
        ac2 = items[0]
        left2 = items[1:] if len(items)>0 else "null"
        ac2_dict[ac2].append(i)
        if ac2 not in acs2:
            acs2.append(ac2)
        for j in items:
            if j[:3] == "IPR":
                if i not in iprs2:
                    ipr = j
                    iprs2_dict[ipr] = i
                    iprs2.append(j)

    #name2 = str(file1[:-4]+"")
    import random
    name = random.choice(acs1)
    try:
        if "/" in  file2:
            namef2 = file2.split("/")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        elif "\\" in file2:
            namef2 = file2.split("\\")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        else:
            namef2 = file2[:-4]
            print namef2
    except:
        print file2
        
    #namef1 = file1.split("/")[-1] 
    namef1 = file1[:-4]
    cf = open(namef1 +"_"+str(len(ac1_dict))+"_"+str(namef2)+"_merge_iprs.txt","w")
    #ccf = open(namef1 +"_"+str(len(ac1_dict))+"_"+str(namef2)+"_merge_iprs2.txt","w")
    #df = open(str(len(ac1_dict))+"_match2_"+namef2+".txt","w")

    cf.write(fields2+"\t"+fields+"\n")
    #cf.write(fields2+"\n")
    
    #df.write(fields+"\n")

    count = 0
    pasted = []
    count2 = 0
    pasted2 = []

    
            

    for ipr, pval in iprs2_dict.iteritems():
        #print  iprs1_dict[ipr] , iprs1_dict[ipr]
        
        #new_line = str(iprs2_dict[ipr])+"\t"+str(iprs1_dict[ipr]
            
        for i, j in iprs1_dict.iteritems():
            
            
            #print i, j
            if ipr == i:
                for a in j:
                    #print i, a
                
                    #for x, y in ac1_dict.iteritems()
                    #for b in iprs2_dict[i]:
                    #print a
                    lines1 = str(iprs2_dict[i])
                    llines  = lines1[2:-2]
                    #print llines
                    tab1 = llines.split("\t")
                    #print tab1
                    line1 = "\t".join(tab1)

                    new_line = str(pval)+"\t"+str(a)+"\n"
                    #print new_line

                    #new_line2 = str(iprs2_dict[i])+"\t"+str(a)
                    #new_line3 = str(iprs2_dict[i])+"\t"+str(a)

                    #print new_line2
                    #print new_line3

                    if new_line not in pasted:
                        pasted.append(new_line)
                        #print i, j
                        cf.write(new_line)
                       
                        #cf.write(new_line+"\n")
                        count +=1
            
                        
    #for i, j in iprs1_dict.iteritems():
    #    if i not in iprs2_dict.keys():
            #for a in ac1_dict[i]:
    #        for a in j:
    #            #print a
     #           new_line = a
#
     #           if new_line not in pasted:
     #               pasted.append(new_line)
                    #print i, j
      #              cf.write(new_line+"\n")
     #               count +=1
            
        

    
    #count2 = 0
    #for i, j in ac2_dict.iteritems():
     #   if i in ac1_dict.keys():
      #      df.write(str(ac2_dict[i])+"\n")
       #     count2 +=1
    print "lines : ", count# , count2
    print "acs1: ", len(ac1_dict)
    print "acs2: ", len(ac2_dict)
    print "TOTAL Lines to New File : ", count
    
    print cf.name, " saved!"
    #print df, " save!"

    theInputFiles.append(cf.name)
    #theInputFiles.append(ccf.name)
    #theInputFiles.append(df.name)



#merge_ac_with_ipr_ipr("20151013_meegeds/tax5671_bmq_ac_ipr_go_k.txt","20151013_meegeds/top_signi.txt")

#merge2files_by_ac("uniprot-taxonomy%3A5671.list", "output/5671_bs_uniprot.txt")
#merge2files_by_ac("tests/2pep_260.txt","tests/2peps.txt")
#output\\A4HZ66_3990_bs_uniprot.txt

#match2files_by_ac("2pep_246_down.txt", "output/tax5671_bmq_ac_ipr_go.txt")

#############################################################
########################################################################
###  AUX MERGE UNI IDS WITH TAX ID FILE


def match_uni_ids_from_list(uni_ids, gene_ids, thefilename):

    global theInputFiles

    af = open(thefilename, "r")

    if len(uni_ids) > 0:
        name = random.choice(uni_ids)
        bf = open(str(thefilename[:-4])+"_filter_"+str(len(uni_ids))+"_"+str(name)+".txt", "w")


        bf.write(af.readline())
        for line in af.readlines():
            line = line.strip()
            items = line.split("\t")
            ac = items[0]
            #print ac
            if ac in uni_ids:
                bf.write(line+"\n")
                #print line

        bf.close()
        af.close()
        print bf, " saved!"
        theInputFiles.append(bf.name)
        print theInputFiles
        
    elif len(gene_ids) > 0:
        name = random.choice(gene_ids)
        bf = open(str(thefilename[:-4])+"_filter_"+str(len(gene_ids))+"_"+str(name)+".txt", "w")


        bf.write(af.readline())
        for line in af.readlines():
            line = line.strip()
            items = line.split("\t")
            
            #print ac
            for i in gene_ids:
                for j in items:
                    item = j.split(" ")
                    if  i in item:
                        bf.write(line+"\n")
                    #print line

        bf.close()
        af.close()
        print bf, " saved!"
        theInputFiles.append(bf.name)
        print theInputFiles
        
    else:
        print "No ID was found in file: ", af.name
        return
        

    #return bf
            
        


def from_ipr_get_ac_infile(ipr_ids, thefilename):

    global theInputFiles

    af = open(thefilename, "r")

    name = random.choice(ipr_ids)
    bf = open(str(thefilename[:-4])+"_filter_"+str(len(ipr_ids))+"_"+str(name)+"_IPR2ACs.txt", "w")
    uniqlines = []
    bf.write(af.readline())
    for line in af.readlines():
        line = line.strip()
        items = line.split("\t")
        #ac = items[0]
        #ipr = items[1] if len(line) >1 else "null"
        
        for i in items:
            if i[:3] == "IPR":
                #print i
                i = i.split("; ")
                #print i
                for a in i:
                    #print a
                    if a in ipr_ids:
                        if line not in uniqlines:
                            bf.write(line+"\n")
                            uniqlines.append(line)#print line
                    
                    

    bf.close()
    af.close()
    print bf, " saved!"
    theInputFiles.append(bf.name)

    #return bf




def from_go_get_ac_infile(go_ids, thefilename):

    global theInputFiles

    af = open(thefilename, "r")

    name = random.choice(go_ids)
    bf = open(str(thefilename[:-4])+"_filter"+str(len(go_ids))+"_"+str(name)+"_GO2ACs.txt", "w")

    bf.write(af.readline())
    print af.readline()
    for line in af.readlines():
        line = line.strip()
        items = line.split("\t")
        ac = items[0]
        go = items[7] if len(line) >7 else "null"
        if go == "GO:":
            #print ac
            if go in go_ids:
                bf.write(line+"\n")
                #print line
        else:
            for i in items:
                if i[:3] == "GO:": # pos[4]
                    i = i.split("; ")
                    for a in i:
                        if a in go_ids:
                            bf.write(line+"\n")
                            #print line
                    
                    

    bf.close()
    af.close()
    print bf, " saved!"
    theInputFiles.append(bf.name)



#####
#word = "kinase"
#filename = "output/all_back.txt"
# regular Patterns :
# http://www.dotnetperls.com/re

def find_word_in_line(filename, words):

    global theInputFiles

    print "\nFind word in File : "
    #word = raw_input("Enter word : ")

    values = words # ["cat100", "---200", "xxxyyy", "jjj", "box4000", "tent500"]
    print words
    
    # Find word in line and write in newfile
    cnt = Counter()
    new_file = open(str(filename[:-4])+"_find_"+str(words[0])+".txt","w")

    af = open(filename,"r")
    header = af.readline()
    new_file.write(header)
    count = 0
    for i in words[:3]: # only 3 words
        words.append(i.lower())
        words.append(i.upper())

    lines = []
    
    with open(filename,"r") as f1:
        for i in f1.readlines():
            line = i.strip()
            items = line.split("\t")
            for j in items:
                #print j
                #j = j.split(" ")
                #for k in j: 
                for word in words:
                    if word in j:
                        if line not in lines:
                            #print j, line
                            new_file.write(line+"\n")
                            count += 1
                            lines.append(line)
                #continue

            
    print new_file.name, " saved!\nComplete!"
    
    theInputFiles.append(new_file.name)



"""
import re
search="please help me out"
fullstring="please help me out so that I could solve this"
s = re.search(search,fullstring)
print(s.group())
"""
#find_most_common() 
#find_words()           
#find_word_in_line(filename,word)

import itertools

def slicefile(filename, start, end):
    lines = open(filename)
    return itertools.islice(lines, start, end)

def cut_my_file(filename, start,end):

    out = open("output\\new_file.txt", "w")
    for line in slicefile(filename, start, end):
        out.write(line)
    out.close()
    print out.name, " saved !"
            
    
#cut_my_file(filename, 10,15)


#####
#filename = "output/all_back.txt"


def split_file():

    with open(filename) as f:
        contents1, sentinel, contents2 = f.read().partition("kinase\n")
    with open(str(filename[:-4]+"_p1.txt"), "w") as f:
        f.write(contents1)
    with open(str(filename[:-4]+"_p2.txt"), "w") as f:
        f.write(contents2)

        
#####

    
# random names
#import uuid
#filename = str(uuid.uuid4())

#import datetime
#basename = "mylogfile"
#suffix = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
#filename = "_".join([basename, suffix]) # e.g. 'mylogfile_120508_171442'

# MAKE FILE DIR
# try/except
#filename = "/foo/bar/baz.txt"
#if not os.path.exists(os.path.dirname(filename)):
#    os.makedirs(os.path.dirname(filename))
#with open(filename, "w") as f:
#    f.write("FOOBAR")



###########################################33    MATCH    ##################################3

from collections import defaultdict

def match2files_by_ac(file1, file2):

    global theInputFiles

    uni_ids = []
    

    af = open(file1, "r") # choice 1
    bf = open(file2,"r") # background


    ac1_dict = defaultdict(list)
    spacer = 77*"#"
    print spacer, "\n\n\tMatching AC IDS \n\t" ,file1 ," VS. ", file2
    acs1 = []
    for i in af.readlines():
        i = i.strip()
        items = i.split("\t")
        ac1 = items[0]
        
        left = items[1:] if len(items)>0 else "null"
        ac1_dict[ac1].append(i)
        if ac1 not in acs1:
            acs1.append(ac1)


    ac2_dict = defaultdict(list)
    acs2 = []
    fields =  bf.readline()
    print fields
    
    for i in bf.readlines()[1:]:
        i = i.strip()
        items = i.split("\t")
        ac2 = items[0]
        left2 = items[1:] if len(items)>0 else "null"
        ac2_dict[ac2].append(i)
        if ac2 not in acs2:
            acs2.append(ac2)

    #name2 = str(file1[:-4]+"")
    name = random.choice(acs1)
    try:
        if "/" in  file2:
            namef2 = file2.split("/")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        elif "\\" in file2:
            namef2 = file2.split("\\")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        else:
            namef2 = file2[:-4]
            print namef2
    except:
        print file2

    #namef1 = file1.split("/")[-1] 
    namef1 = file1[:-4]
    cf = open(str(namef1)+"_filter_"+str(len(ac1_dict))+"_"+str(namef2)+".txt","w")
    #df = open(str(len(ac1_dict))+"_match2_"+namef2+".txt","w")

    cf.write(fields)
    #df.write(fields+"\n")

    count = 0
    pasted = []
    for i, j in ac2_dict.iteritems():
        if i in ac1_dict.keys():
            for a in j:
                if a not in pasted:
                    pasted.append(a)
                    #print i, j
                    cf.write(str(a)+"\n")
                    count +=1

    #count2 = 0
    #for i, j in ac2_dict.iteritems():
     #   if i in ac1_dict.keys():
      #      df.write(str(ac2_dict[i])+"\n")
       #     count2 +=1
    print "lines : ", count# , count2
    print "acs1: ", len(ac1_dict)
    print "acs2: ", len(ac2_dict)
    print cf, " save!"
    #print df, " save!"

    theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)

    #thefilename = from_tax_id_get_bs_uniprot(new_tax_id)
    #match_uni_ids_from_list(uni_ids, thefilename.name)

###############################################################################
    
def merge_inline(file1, file2):

    global theInputFiles
    from collections import defaultdict

    info = "Select bigger file ( ComboBox 1) for yield matching. File 1 in Column 1 ."


    uni_ids = []
    
    af = open(file1, "r") # choice 1
    bf = open(file2,"r") # background

    ac1_dict = defaultdict(list)
    spacer = 77*"#"
    print spacer, "\n", info, "\n\n\tMatching AC IDS \n\t" ,file1 ," VS. ", file2
    
    acs1 = []
    fields =  af.readline()
    fields = fields.strip()
    for i in af.readlines():
        i = i.strip()
        items = i.split("\t")
        ac1 = items[0]
        
        left = items[1:] if len(items)>0 else "null"
        ac1_dict[ac1].append(i)
        if ac1 not in acs1:
            acs1.append(ac1)


    ac2_dict = defaultdict(list)
    acs2 = []
    fields2 =  bf.readline()
    fields2 = fields2.strip()
    
    print fields, "\t", fields2
    
    for i in bf.readlines()[1:]:
        i = i.strip()
        items = i.split("\t")
        ac2 = items[0]
        left2 = items[1:] if len(items)>0 else "null"
        ac2_dict[ac2].append(i)
        if ac2 not in acs2:
            acs2.append(ac2)

    #name2 = str(file1[:-4]+"")
    import random
    name = random.choice(acs1)
    try:
        if "/" in  file2:
            namef2 = file2.split("/")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        elif "\\" in file2:
            namef2 = file2.split("\\")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        else:
            namef2 = file2[:-4]
            print namef2
    except:
        print file2
        
    #namef1 = file1.split("/")[-1] 
    namef1 = file1[:-4]
    cf = open(namef1 +"_"+str(len(ac1_dict))+"_"+str(namef2)+"_merged.txt","w")
    #df = open(str(len(ac1_dict))+"_match2_"+namef2+".txt","w")

    cf.write(fields+"\t"+fields2+"\n")
    #cf.write(fields2+"\n")
    
    #df.write(fields+"\n")

    count = 0
    pasted = []
    for i, j in ac1_dict.iteritems():
        if i in ac2_dict.keys():
            #for x, y in ac1_dict.iteritems()
            for b in ac2_dict[i]:
                for a in j:
                    #print a
                    new_line = a+"\t"+b
                    
                    if new_line not in pasted:
                        pasted.append(new_line)
                        #print i, j
                        cf.write(new_line+"\n")
                        count +=1
    for i, j in ac1_dict.iteritems():
        if i not in ac2_dict.keys():
            #for a in ac1_dict[i]:
            for a in j:
                #print a
                new_line = a

                if new_line not in pasted:
                    pasted.append(new_line)
                    #print i, j
                    cf.write(new_line+"\n")
                    count +=1
            
        

    
    #count2 = 0
    #for i, j in ac2_dict.iteritems():
     #   if i in ac1_dict.keys():
      #      df.write(str(ac2_dict[i])+"\n")
       #     count2 +=1
    print "lines : ", count# , count2
    print "acs1: ", len(ac1_dict)
    print "acs2: ", len(ac2_dict)
    print "TOTAL Lines to New File : ", count
    
    print cf.name, " saved!"
    #print df, " save!"

    theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)





#merge2files_by_ac("uniprot-taxonomy%3A5671.list", "output/5671_bs_uniprot.txt")
#merge2files_by_ac("tests/2pep_260.txt","tests/2peps.txt")
#output\\A4HZ66_3990_bs_uniprot.txt

#match2files_by_ac("2pep_246_down.txt", "output/tax5671_bmq_ac_ipr_go.txt")

#############################################################

def match2files_by_ac_only(file1, file2):

    global theInputFiles

    uni_ids = []
    

    af = open(file1, "r") # choice 1
    bf =open(file2,"r") # background


    ac1_dict = defaultdict(list)
    spacer = 77*"#"
    print spacer, "\n\n\tMatching by AC IDS \n\t" ,file1 ," + ", file2
    acs1 = []
    for i in af.readlines():
        i = i.strip()
        items = i.split("\t")
        ac1 = items[0]
        
        left = items[1:] if len(items)>0 else "null"
        ac1_dict[ac1].append(i)
        if ac1 not in acs1:
            acs1.append(ac1)


    ac2_dict = defaultdict(list)
    acs2 = []
    fields =  bf.readline()
    print fields
    
    for i in bf.readlines():
        i = i.strip()
        items = i.split("\t")
        ac2 = items[0]
        left2 = items[1:] if len(items)>0 else "null"
        ac2_dict[ac2].append(i)
        if ac2 not in acs2:
            acs2.append(ac2)

    #name2 = str(file1[:-4]+"")
    name = random.choice(acs1)
    try:
        if "/" in  file2:
            namef2 = file2.split("/")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        elif "\\" in file2:
            namef2 = file2.split("\\")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        else:
            namef2 = file2[:-4]
            print namef2
    except:
        print file2
        
    #namef1 = file1.split("/")[-1] 
    namef1 = file1[:-4]
    cf = open(str(namef1)+"_"+str(len(ac1_dict))+"_"+str(namef2)+"_match.txt","w")
    #df = open(str(len(ac1_dict))+"_match2_"+namef2+".txt","w")

    cf.write(fields+"\n")
    #df.write(fields+"\n")

    count = 0
    pasted = []
    nomatch = []
    for i, j in ac1_dict.iteritems():
        for k,v in ac2_dict.iteritems():
            if i == k: # ONLY WHEN MATCH BOTH FILES
                for a in j:
                    if a not in pasted:
                        pasted.append(a)
                        #print i, j
                        cf.write(str(a)+"\n")
                        count +=1
                for b in v:
                    if b not in pasted:
                        pasted.append(b)
                        cf.write(str(b)+"\n")
                        count +=1

    for i,j in ac1_dict.iteritems():
        if i not in ac2_dict.keys():
            nomatch.append(i)
            for a in j:
                cf.write(str(a)+"\n")
        
                
            
    #count2 = 0
    #for i, j in ac2_dict.iteritems():
     #   if i in ac1_dict.keys():
      #      df.write(str(ac2_dict[i])+"\n")
       #     count2 +=1
       
   
        
        
    print "lines : ", count# , count2
    print "acs1: ", len(ac1_dict)
    print "acs2: ", len(ac2_dict)
    print cf, " save!"
    #print df, " save!"

    theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)


##########################################################3

def merge2files_by_ac_all(file1, file2):

    global theInputFiles

    uni_ids = []
    

    af = open(file1, "r") # choice 1
    bf =open(file2,"r") # background


    ac1_dict = defaultdict(list)
    spacer = 77*"#"
    print spacer, "\n\n\tMerging by AC IDS \n\t" ,file1 ," + ", file2
    acs1 = []
    for i in af.readlines():
        i = i.strip()
        items = i.split("\t")
        ac1 = items[0]
        
        left = items[1:] if len(items)>0 else "null"
        ac1_dict[ac1].append(i)
        if ac1 not in acs1:
            acs1.append(ac1)


    ac2_dict = defaultdict(list)
    acs2 = []
    fields =  bf.readline()
    print fields
    
    for i in bf.readlines():
        i = i.strip()
        items = i.split("\t")
        ac2 = items[0]
        left2 = items[1:] if len(items)>0 else "null"
        ac2_dict[ac2].append(i)
        if ac2 not in acs2:
            acs2.append(ac2)

    #name2 = str(file1[:-4]+"")
    name = random.choice(acs1)
    try:
        if "/" in  file2:
            namef2 = file2.split("/")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        elif "\\" in file2:
            namef2 = file2.split("\\")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        else:
            namef2 = file2[:-4]
            print namef2
    except:
        print file2
        
    #namef1 = file1.split("/")[-1] 
    namef1 = file1[:-4]
    cf = open(str(namef1)+"_"+str(len(ac1_dict))+"_"+str(namef2)+"_merged.txt","w")
    #df = open(str(len(ac1_dict))+"_match2_"+namef2+".txt","w")

    cf.write(fields+"\n")
    #df.write(fields+"\n")

    count = 0
    pasted = []
    for i, j in ac1_dict.iteritems():
        for k,v in ac2_dict.iteritems():
            if i == k: # ONLY WHEN MATCH BOTH FILES
                for a in j:
                    if a not in pasted:
                        pasted.append(a)
                        #print i, j
                        cf.write(str(a)+"\n")
                        count +=1
                for b in v:
                    if b not in pasted:
                        pasted.append(b)
                        cf.write(str(b)+"\n")
                        count +=1
            
 
    for i, j in ac1_dict.iteritems():
        if i not in ac2_dict.keys():
            for a in j:
                if a not in pasted:
                    pasted.append(a)
                    #print i, j
                    cf.write(str(a)+"\n")
                    count +=1

    for i, j in ac2_dict.iteritems():
        if i not in ac1_dict.keys():
            for a in j:
                if a not in pasted:
                    pasted.append(a)
                    #print i, j
                    cf.write(str(a)+"\n")
                    count +=1
    
            
                    

    #count2 = 0
    #for i, j in ac2_dict.iteritems():
     #   if i in ac1_dict.keys():
      #      df.write(str(ac2_dict[i])+"\n")
       #     count2 +=1
    print "lines : ", count# , count2
    print "acs1: ", len(ac1_dict)
    print "acs2: ", len(ac2_dict)
    print cf, " save!"
    #print df, " save!"

    theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)

    #thefilename = from_tax_id_get_bs_uniprot(new_tax_id)
    #match_uni_ids_from_list(uni_ids, thefilename.name)


#merge2files_by_ac("uniprot-taxonomy%3A5671.list", "output/5671_bs_uniprot.txt")
#merge2files_by_ac("tests/2pep_260.txt","tests/2peps.txt")
#output\\A4HZ66_3990_bs_uniprot.txt

#merge2files_by_ac_all("output/dynein.txt", "tes15.txt")


def join_lines(file1, file2):

    global theInputFiles

    uni_ids = []
    

    af = open(file1, "r+") # choice 1
    bf =open(file2, "r+") # background

    count = 0
    pasted = []
    nomatch = []
    
    fields =  af.readline()
    print fields

    #name2 = str(file1[:-4]+"")
    
    try:
        if "/" in  file2:
            namef2 = file2.split("/")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        elif "\\" in file2:
            namef2 = file2.split("\\")
            namef2 = namef2[-1]
            namef2 = namef2[:-4]
        else:
            namef2 = file2[:-4]
            #print namef2
    except:
        print file2
        
    #namef1 = file1.split("/")[-1] 
    namef1 = file1[:-4]
    cf = open(str(namef1)+"_"+str(namef2)+"_join.txt","w")
    #df = open(str(len(ac1_dict))+"_match2_"+namef2+".txt","w")

    cf.write(fields+"\n")
    #df.write(fields+"\n")

    print " Will remove repeated lines. . . "


    ac1_dict = defaultdict(list)
    spacer = 77*"#"
    print spacer, "\n\n\tJoining Lines \n\t" ,file1 ," + ", file2

    acs1 = []
    for i in af.readlines():
        i = i.strip()
        if i not in pasted:
            pasted.append(i)
            cf.write(i+"\n")
            count +=1
            
        items = i.split("\t")
        ac1 = items[0]
        
        left = items[1:] if len(items)>0 else "null"
        ac1_dict[ac1].append(i)
        if ac1 not in acs1:
            acs1.append(ac1)


    ac2_dict = defaultdict(list)
    acs2 = []
    
    
    for i in bf.readlines():
        i = i.strip()

        if i not in pasted:
            pasted.append(i)
            cf.write(i+"\n")
            count += 1
            
        items = i.split("\t")
        ac2 = items[0]
        
        left2 = items[1:] if len(items)>0 else "null"
        ac2_dict[ac2].append(i)
        if ac2 not in acs2:
            acs2.append(ac2)

    
   
    print "acs1: ", len(ac1_dict)
    print "acs2: ", len(ac2_dict)
    print "lines : ", count# , count2
    
    print cf, " save!"
    #print df, " save!"

    theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)


#join_lines("output/dynein.txt", "tes15.txt")


########################### AUX ##################################
#### SAVE UNIQUE
####################################################################
# All files have accession number (ac) at left position in the table
# All files have taxonomy number (tax) at right position in the table

def save_input_file_to_append_by_ac(theInput1, theInput2, theOutput):
    aFile = open(theInput1, "rt")
    cFile = open(theOutput, "wt")
    bFile = open(theInput2, "rt")
    entries = []
    entries_dict = {}
    count = 0
    for aLine in bFile:
            aProtein_entry2 = aLine.strip().split('\t')
            entry = aProtein_entry2[0]
            entries_dict[aProtein_entry2[0]] = aLine # Assuming ac is in first position
            entries.append(aLine)
    for aLine in aFile:
            aProtein_entry1 = aLine.strip().split('\t')
            if aProtein_entry1[0] in entries_dict.keys():
                    print aProtein_entry1[0], aLine
                    #aProtein_entry1 = aLine.strip().split('\t')
                    #print aProtein_entry[0]
                    cFile.write(aLine)
                    count  += 1
                    
    print "Accessions matched ", count, "times"
    aFile.close()
    bFile.close()
    cFile.close()

"""A4IE50	A4IE50_LEIIN	unreviewed	Uncharacterized protein	LINJ_36_4440		Predicted	1123	1	11	11	18	3,37E+11	100.000	2	0.0	122.7	5.52
"""
#theInput1 = "data1.txt"
#theInput2 = "data2.txt"
#theOutput = raw_input("Save file as : ")
# save_input_file_to_append_by_ac(theInput1, theInput2,theOutput )
########################################################################################33


def filter_entries(theFilename):

    global theInputFiles

    af = open(theFilename,"rt")
    cf = open(str(theFilename[:-4])+"_filtered.txt","wt")
    entries_list = []
    count = 0
    for i in af.readlines():
        count +=1
        if i not in entries_list:
            entries_list.append(i)
            cf.write(i)

    print "\nFile 1 : total  entries = ", count
    print "File 2 : total entries = ", len(entries_list)

    print cf, "save!"
    af.close()
    cf.close()
    theInputFiles.append(cf.name)




###
def merge_same_line(theFilename1,theFilename2):
    
    from collections import defaultdict
    entries1 = []
    entries2 = []
    entries3 = []
    entries_dict1 = defaultdict(list)
    
    af = open(theFilename1,"r")
    bf = open(theFilename2, "w")
    acs = []
    alines = []
    new_line = ""
    for line in af.readlines():
        alines.append(line.strip())
        ac1 = line.split("\t")[0]
        if ac1 not in acs:
            acs.append(ac1)
            #print ac1

    #print "hello"
    for i in acs:
        for line in alines:
            ac = line.split("\t")[0]
            if i == ac:
                #print i
                line = line+"\t"
                #print line
                new_line = new_line + line
            else:
                pass
        #print new_line, "\n\n"
        
        #for i in all_lines:
        bf.write(new_line+"\n")
        new_line = ""

    af.close()
    bf.close()
    print bf.name, "saved!"
    print "complete!"
                
            

#merge_same_line("260_prot_merged_filtered.txt", "260_prot_merged_filtered_line.txt")



# I HAVE 2 FILES TO READ AND ONE I CANNOT CHANGE THE ORDER, MAX LIGHT, FOR EXAMPLE
# ON THE OTHER HAND I CAN SORT THE OTHER LIST
# I AM GOIN TO ASSUME YOUR LIST HAS A SPECIFIC ORDER.

# MERGE BY AC ID

def save_input_file_to_append_by_ac_ids(theInput1, theInput2,theInput3,theInput4,theInput5,theInput6,theInput7, theOutput):
    letters= ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","X","Y","Z"]
    numbers = ["1","2","3","4","5","6","7","8","9","0"]
    if theInput1 is not "theInputList1":
        aFile = open(str(theInput1), "rt")
    if theInput2 is not "theInputList2":
        bFile = open(str(theInput2), "rt")
    #if theInput4 is not None:
     #   eFile = open(str(theInput4), "rt") 
    #if theInput5 is not None:
    #    fFile = open(str(theInput5), "rt") 
    #if theInput6 is not None:
     #   gFile = open(str(theInput6), "rt") 
    cFile = open("output/"+str(theOutput), "wt")
    from collections import defaultdict
    entries1 = []
    entries2 = []
    entries3 = []
    entries_dict1 = defaultdict(list)
    entries_dict2 = defaultdict(list)
    entries_dict3 = defaultdict(list)
    entries_dict4 = defaultdict(list)
    entries_dict5 = defaultdict(list)

    count = 0
    IDs = 0
    theUniIDs = []
    theTaxIDs = []
    lines1 = 0
    lines2 = 0
    lines3 = 0
    for aLine in aFile:
            aProtein_entry1 = aLine.strip().split('\t')
            if len(aProtein_entry1) == 1:
                line = "null"
            else:
                line = aProtein_entry1[0]
                spaces1 = len(aProtein_entry1)
                if line[0] in letters and line[1] in numbers and line[-1] in numbers:                    
                    lines1 += 1
                    entries_dict1[line].append(aLine) # Assuming ac is in first position
                    entries1.append(aLine)
                    # print filter(lambda x: x > 0, entries)
    for aLine in bFile:
            #print aLine
            aProtein_entry2 = aLine.strip().split('\t')
            if len(aProtein_entry2) == 1:
                line = "null"
            else:
                line = aProtein_entry2[0] 
                spaces2 = len(aProtein_entry2)
                if line[0] in letters and line[1] in numbers and line[-1] in numbers:
                    lines2 += 1
                    # data_dict[regNumber].append(details)
                    entries_dict2[line].append(aLine) # Assuming ac is in first position
                    entries2.append(aLine)
                    # print filter(lambda x: x > 0, entries

    if theInput3 is not "theInputList3":
        dFile = open(str(theInput3), "rt") # if len(theInput1) > 0 : print "loading " # "null": print "Load Input 1 "
        for aLine in dFile.readlines():
            #print aLine
            aProtein_entry3 = aLine.strip().split('\t')
            if len(aProtein_entry3) == 1:
                line = "null"
            else:
                line = aProtein_entry3[0]
                spaces3 = len(aProtein_entry3)
                if line[0] in letters and line[1] in numbers and line[-1] in numbers:
                    lines3 += 1
                    # data_dict[regNumber].append(details)
                    entries_dict3[line].append(aLine) # Assuming ac is in first position
                    entries3.append(aLine)
                    # print filter(lambda x: x > 0, entries
        dFile.close()
        
    count_yes = 0
    count_no = 0
    for i in entries_dict1.keys():
        #print entries_dict1[i]
        if i in entries_dict2.keys() and i in entries_dict3.keys():
            #print i
            merged = entries_dict1[i] + entries_dict2[i] + entries_dict3[i] # append more dicts from input
            #merged = "".join(merged)
            #merged = merged.replace("\n","\t")
            #print merged
            #cFile.write(merged+"\n")
            for i in merged:
                i.replace("\n","\t")
                cFile.write(i)
            IDs += 1
                
        elif i in entries_dict2.keys():
            #print i
            merged = entries_dict1[i] + entries_dict2[i]
            #merged = "".join(merged)
            #merged = merged.replace("\n","\t")
            #print merged
            #cFile.write(merged+"\n")
            for i in merged:
                i.replace("\n","\t")
                cFile.write(i)
            IDs += 1
                
        elif i in entries_dict3.keys():
            #print i
            merged = entries_dict1[i] + entries_dict3[i]
            #merged = "".join(merged)
            #merged = merged.replace("\n","\t")
            #print merged
            #cFile.write(merged+"\n")
            for i in merged:
                i.replace("\n","\t")
                cFile.write(i)
            IDs += 1
                
        else:
            count_no += 1
            #print i, "has no match"
            merged = entries_dict1[i]
            #merged = "".join(merged)
            #merged = merged.replace("\n","\t")
            #print merged
            #cFile.write(merged+"\n")
            for i in merged:
                i.replace("\n","\t") #to tryng extend listi
                cFile.write(i)
            IDs += 1




    print "\nfile1 - total entries and accessions : ", len(entries1), len(entries_dict1)
    print "file2 - total entries and accesions : ", len(entries2), len(entries_dict2)
    print "file3 - total entries and accesions : ", len(entries3), len(entries_dict3)
    print "TOTAL IDs appended : ", IDs
    print "TOTAL IDS not appended : ", count_no
    
    aFile.close()
    bFile.close()
    
    
    cFile.close()
    print cFile,"saved!"
    
    filter_entries(theOutput)
    #theOut = raw_input("Save file as : ")
    cfile_short = cFile.name[:-4]
    #print cfile_short
    merge_same_line(cFile.name, str(cfile_short)+"_line.txt")
    
    print "(complete!)"


   



"""5671	ATAT_LEIIN	A4I1F7	IPR016181; IPR007965; 	Alpha-tubulin N-acetyltransferase (Alpha-TAT) (TAT) (EC 2.3.1.108) (Acetyltransferase mec-17 homolog)	LinJ25.1190 LinJ_25_1190	GO:0071929; GO:0070507; GO:0019799	alpha-tubulin acetylation; regulation of microtubule cytoskeleton organization; tubulin N-acetyltransferase activity		Inferred from homology	Leishmania infantum
A4I1F7	ATAT_LEIIN	reviewed	Alpha-tubulin N-acetyltransferase (Alpha-TAT) (TAT) (EC 2.3.1.108) (Acetyltransferase mec-17 homolog)	LinJ25.1190 LinJ_25_1190	GO:0071929; GO:0070507; GO:0019799	Inferred from homology	246	1	2	2	2	6.568E7	5.608	1		27.1	9.61
"""

#theInputFile1  =  "2pep_10.txt"
#theInputFile2  =  "output/E9AHJ2_bmq_ac_ipr_go_db.txt"
#theInputFile3  =  "kegg_data/E9AHJ2_get_kegg_paths.txt"
#theInputFile4 = "data1.txt" # 5 3
#theInputFile5 = "data2.txt" # 244 244 ho to extend=?
#theInputFile6 = "data0.txt" # 11 3
#theInputFile7 = "5671_ac_test.txt" # 14

#theOutput = str(raw_input("Save file as : ")+"_merged.txt")
#save_input_file_to_append_by_ac_ids(theInputFile1, theInputFile2,theInputFile3,theInputFile4 ,theInputFile5 ,theInputFile6 ,theInputFile7 ,theOutput )

###
def remove_duplicates():
    t = ['a', 'b', 'c', 'd']
    t2 = ['a', 'c', 'd']
    for t in t2:
        t.append(t.remove())
    return t


def remove_exponent(value):
    """
       >>>(Decimal('5E+3'))
       Decimal('5000.00000000')
    """
    decimal_places = 8
    max_digits = 16

    if isinstance(value, decimal.Decimal):
        context = decimal.getcontext().copy()
        context.prec = max_digits
        return "{0:f}".format(value.quantize(decimal.Decimal(".1") ** decimal_places, context=context))
    else:
        return "%.*f" % (decimal_places, value)



###############################################################################################3


def fileinput(files):
    for f in files:
        with open('r') as fin:
            for line in fin:
                yield line

#bar_plot2labels(all_files)      
# bar_plot2labels(af,bf,cf,df,ef,kf) 


###### ################      CONVERT TXT tO EXCEL     ############################

# import os
# os.rename("C:\Users\JohnDoe\Downloads\Report.txt", "C:\Users\JohnDoe\Downloads\Report.cvs")

import Tkinter
import tkFileDialog
from tkMessageBox import *

def convert_txt2excel(textfile):

    # https://www.daniweb.com/software-development/python/threads/467101/convert-text-files-to-excel-files-

    print "################   TEXT TO EXCEL CONVERSION   #################\n"

    from os import listdir
    from os.path import isfile, join
    import xlwt
    import xlrd

    global theInputFiles

    # mypath should be the complete path for the directory containing the input text files
    # mypath = 'C:\\Temp'
    #mypath = "C:\Users\mendes\Dropbox\Thesis\aeferreira\leishdomains\output"

    #ff = tkFileDialog.asksaveasfile(mode='w', defaultextension=".xls") ERROR
    #ff = tkFileDialog.asksaveasfilename( defaultextension=".xls")
    #if ff is None: # asksaveasfile return `None` if dialog closed with "cancel".
    #return
    #name = os.path.basename(ff)

    #textfiles = [ join(mypath,f) for f in listdir(mypath) if isfile(join(mypath,f)) and '.txt' in  f]
    #for textfile in textfiles:
    
    f = open(textfile, 'r+')
    row_list = []
    for row in f:
        row_list.append(row.split('\t'))
    column_list = zip(*row_list)
    workbook = xlwt.Workbook()
    i = 0
    sheetname = "Sheet%s" % i
    #int sheetname
    worksheet = workbook.add_sheet(sheetname)
    
    for column in column_list:
        for item in range(len(column)):
            worksheet.write(item, i, column[item])
        i+=1

    book_name = textfile[:-4] + '.xls'
    #book_name = name

    workbook.save(book_name)
    namepath = os.path.abspath(f.name) #+ book_name
    mypath = namepath[:-4]+   '.xls'
    showinfo(" TXT2EXCEL Conversion ", "File Saved :\n"+mypath)

    print book_name, " saved!"
    theInputFiles.append(book_name)

    print "(complete!)"

#convert_txt2excel("test34.txt")

################################################################################


def convertall_txt2excel(theInputs):

    # https://www.daniweb.com/software-development/python/threads/467101/convert-text-files-to-excel-files-

    print "################   TEXT TO EXCEL CONVERSION  ALL  #################\n"

    from os import listdir
    from os.path import isfile, join
    import xlwt
    import xlrd

    global theInputFiles
    
    # mypath should be the complete path for the directory containing the input text files
    # mypath = 'C:\\Temp'
    #mypath = "C:\Users\mendes\Dropbox\Thesis\aeferreira\leishdomains\output"

    #ff = tkFileDialog.asksaveasfile(title="Save Sesion File as ",mode='w', defaultextension=".xls")

    ff = tkFileDialog.asksaveasfilename( defaultextension=".xls")
    if ff is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return

    names1 = glob('*.xls')
    names2 = glob('*/*.xls')
    
    namepath = os.path.abspath(ff)
    

    namebase = os.path.basename(namepath)
    namedir = os.path.dirname(namepath)
    

    folder = namedir.split("\\")[-1]

    namewdir =folder + "\\" + namebase
    
    #print namepath
    #textfiles = [ join(mypath,f) for f in listdir(mypath) if isfile(join(mypath,f)) and '.txt' in  f]
    #for textfile in textfiles:
        
    row_list = []
    workbook = xlwt.Workbook()
    print theInputs
    
    if ff:
        for j,k in enumerate(theInputs):
            #fname = k.split("\\")[-1]
            print j, k
            f = open(k, 'r')
            sheetname = "Sheet%s" % j
            #print sheetname
            worksheet = workbook.add_sheet(sheetname)
            i = 0

            for row in f:
                row_list.append(row.split('\t'))
                column_list = zip(*row_list)

            for column in column_list:
                for item in range(len(column)):
                    worksheet.write(item, i, column[item])
                i+=1

            thepath = os.path.abspath(f.name)
            worksheet.write(0,i,"File_info")
            worksheet.write(1,i,thepath)
            worksheet.write(2,i,k)
            row_list = []

            #book_name = name
            #book_name = "babe"+".xls"
            workbook.save(namepath)
            showinfo(" Excel Conversion Export ", "File Saved :\n"+namepath)

            print namepath, " saved!\n"
            #print os.path.abspath(workbook)
            if namebase in names1:
                theInputFiles.append(namebase)
            elif namewdir in names2:
                theInputFiles.append(namewdir)

            print "(complete!)"
    else:
        showinfo(" Excel Export All ", "No Excel was selected.")            
        

#convert_txt2excel("tes55.txt")
#theInputFiles = ["20150831_estes\\2peps_3999_filter_207_A4HTV2.txt","20150831_updownsums\\tax5671_bs_uniprot_filter_1088_A4I0E7_counts_ipr.txt"]
   
#convertall_txt2excel(theInputFiles)


##################################   TXT 2 EXCEL FORMAT   3######################

def is_number(s):

    try:
        float(s)
        return True
    except ValueError:
        return False
    
def excelworksheet():

    # mypath should be the complete path for the directory containing the input text files
    mypath = raw_input("Please enter the directory path for the input files: ")

    from os import listdir
    from os.path import isfile, join
    textfiles = [ join(mypath,f) for f in listdir(mypath) if isfile(join(mypath,f)) and '.txt' in  f]


    import xlwt
    import xlrd

    global theInputFiles

    style = xlwt.XFStyle()
    style.num_format_str = '#,###0.00'

    for textfile in textfiles:
        f = open(textfile, 'r+')
        row_list = []

        for row in f:
            row_list.append(row.split('|'))
        column_list = zip(*row_list)
        workbook = xlwt.Workbook()
        worksheet = workbook.add_sheet('Sheet1')
        i = 0
        for column in column_list:
            for item in range(len(column)):
                value = column[item].strip()
                if is_number(value):
                    worksheet.write(item, i, float(value), style=style)
                else:
                    worksheet.write(item, i, value)
            i+=1
        workbook.save(textfile.replace('.txt', '.xls'))

    print "(complete!)"


###############################################################3
def onSaveImage(self):
    ftypes = [('PNG', '*.png'), ('JPEG', '*.jpg'), ('PDF', '*.pdf')]
    imagefname = tkFileDialog.asksaveasfilename(parent=self, filetypes = ftypes)

