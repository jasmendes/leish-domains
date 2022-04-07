# -*- coding: cp1252 -*-

# BIOSERVICES - Result Analysis, on 2-Jun-2014
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

from matplotlib import pyplot as plt
from matplotlib_venn import *#venn3, venn3_circles module venn##""!EQW

import collections
from collections import Counter

import numbers
import decimal

from scipy import stats


from bioservices_functions_tk import *

global theInputFiles


import winsound
    


##########################################################################
#### UNIPROT - BIOMART 
##########################################################################
#INPUT: File_Format<
#'5671
#'LIPA_LEIIN',
#'A4HY57',
#'IPR013785; IPR006638; IPR003698; IPR007197; ',
#'Lipoyl synthase, mitochondrial (EC 2.8.1.8) (Lipoate synthase) (LS) (Lip-syn) (Lipoic acid synthase)',
#'LinJ19.0190 LinJ_19_0350',
#'GO:0051539; GO:0016992; GO:0046872; GO:0005739',
#'4 iron, 4 sulfur cluster binding; lipoate synthase activity; metal ion binding; mitochondrion',
#'Protein modification; protein lipoylation via endogenous pathway; protein N(6)-(lipoyl)lysine from octanoyl-[acyl-carrier-protein]: step 2/2. ',
#'Mitochondrion. ',
#'Inferred from homology',
#'Leishmania infantum'>

"""[0]Organism ID
[1]Entry name
[2]Entry
[3]InterPro
[4]Protein names
[5]Gene names
[6]Gene ontology IDs
[7]Gene ontology (GO)
[8]Pathway
[9]Subcellular location
[10]Protein existence
[11]Organism"""

"""
[0] Entry
[1] Entry Name
[2] Protein names
[3] InterPro
[4] Gene ontology IDs
[5] Gene ontology (GO)
[6] Gene names
[7] Pathway
[8] Subcellular Location
[9] Protein existence
[10] Organism
[11] Organism ID
"""



#Renato: trabalho Est�gio Fase0
"""Entry
Entry name
Status
Protein names
Gene names
Gene ontology IDs
Protein existence
Length"""

def readUniProtTab (filename, n = 0, with_fields = False):
    f = open (filename)
    first = f.readline()
    first = first.strip()
    fields = first.split('\t')
    protdata = []

    c = 0
    for i in f:
        i = i.strip()
        if len(i) == 0: continue
        data = i.split('\t')
        dprot = {}
        for k, d in zip(fields, data):
            #print k
            dprot[k] = d
        protdata.append(dprot)
        #print dprot
        if n > 0:
            c +=1
            if c == n: break
    f.close()
    if not with_fields:
        return protdata
    else:
        return protdata, fields



####################################################################### 1 april 2015


def readUniProtTab_sum (filename, n = 0, with_fields = False):
    f = open (filename)

    fields = []
    dict_values = {}
    while len(fields) <=1:
        first = f.readline()
        first = first.strip()
        print (first)
        fields = first.split("\t")
        if len(fields) >1:
            for i, value in enumerate(fields):
                print ("col%s =%s " % (i ,value))
                dict_values[i] = value
        break

    #print 
    protdata = []

    one = []
    two = []
    three = []
    four = []
    five = []
    six = []
    seven = []
    eight = []
    nine = []
    ten = []
    eleven = []
    twelve = []
    thirteen = []
    fourteen = []
    fifteen = []
    sixteen = []
    seventeen = []
    eighteen = []
    nineteen = []
    twenty = []
    
    dprot= {}

    c = 0
    for i in f:
        i = i.strip()
        
        if len(i) == 0: continue

        else:
            data = i.split('\t') 
            one.append(data[0]) if len(data) >= 1 else 'null'
            two.append(data[1]) if len(data) >= 2 else 'null'
            three.append(data[2]) if len(data) >= 3 else 'null'
            four.append(data[3]) if len(data) >= 4 else 'null'
            five.append(data[4]) if len(data) >= 5 else 'null'
            six.append(data[5]) if len(data) >= 6 else 'null'
            seven.append(data[6]) if len(data) >= 7 else 'null'
            eight.append(data[7]) if len(data) >= 8 else 'null'
            nine.append(data[8]) if len(data) >= 9 else 'null'
            ten.append(data[9]) if len(data) >= 10 else 'null'
            eleven.append(data[10]) if len(data) >= 11 else 'null'
            twelve.append(data[11]) if len(data) >= 12 else 'null'
            thirteen.append(data[12]) if len(data) >= 13 else 'null'
            fourteen.append(data[13]) if len(data) >= 14 else 'null'
            fifteen.append(data[14]) if len(data) >= 15 else 'null'
            sixteen.append(data[15]) if len(data) >= 16 else 'null'
            seventeen.append(data[16]) if len(data) >= 17 else 'null'
            eighteen.append(data[17]) if len(data) >= 18 else 'null'
            nineteen.append(data[18]) if len(data) >= 19 else 'null'
            twenty.append(data[19]) if len(data) >= 20 else 'null'
            #ten.append(data[9])
            #eleven.append(data[10])
            #twelve.append(data[11]) 
        
            dprot = {}
            for k, d in zip(fields, data):
                #print k
                dprot[k] = d
            protdata.append(dprot)
        #print dprot
        #for i in data:
         #   print data
        if n > 0:
            c +=1
            if c == n: break
    f.close()
    
    all_lists = []
    if len(one) > 0:
        all_lists.append(one)
    if len(two) > 0:
        all_lists.append(two)
    if len(three) > 0:
        all_lists.append(three)
    if len(four) > 0:
        all_lists.append(four)
    if len(five) > 0:
        all_lists.append(five)
    if len(six) > 0:
        all_lists.append(six)
    if len(seven) > 0:
        all_lists.append(seven)
    if len(eight) > 0:
        all_lists.append(eight)
    if len(nine) > 0:
        all_lists.append(nine)
    if len(ten) > 0:
        all_lists.append(ten)
    if len(eleven) > 0:
        all_lists.append(eleven)
    if len(twelve) > 0:
        all_lists.append(twelve)
    if len(thirteen) > 0:
        all_lists.append(thirteen)
    if len(fourteen) > 0:
        all_lists.append(fourteen)
    if len(fifteen) > 0:
        all_lists.append(fifteen)
    if len(sixteen) > 0:
        all_lists.append(sixteen)
    if len(seventeen) > 0:
        all_lists.append(seventeen)
    if len(eighteen) > 0:
        all_lists.append(eighteen)
    if len(nineteen) > 0:
        all_lists.append(nineteen)
    if len(twenty) > 0:
        all_lists.append(twenty)
   

    print ("\n > > > > > > > > > > > > > > > SUMMARY < < < < < < < < < < < < < < < ")
    list3 = fields
    #list3.append("LAA")
    #list3.append("LEE")
    new_file_name = f.name[:-4]
    jfile = open(str(new_file_name)+"_sum_counts.txt", "wt")
    headerj = "\t".join(fields)
    jfile.write(headerj)
    count = 0
    for idx,item in enumerate(all_lists):  # count in range(0,len(all_lists)):
        a = Counter(item)
        col = list3.pop(0)
        #list3.append(col)
        print ("\nTOTAL %s found : %s "  % (col,len(item))) #len(accs) # 6597
        print ("UNIQUE %s :%s " % (col, len(a)) )
        #cov = (len(a)/len(item))
        div = a.most_common(1)
        for i,j in div:
            #print i, "#",j
            #print type(i) # str
            if type(i) == type(1.0):#s_number(i) == True:
                #isinstance(x, (int, long, float, complex)) == True:
                print ("\n")
            else:
                cov = 1 - (j/float(len(item))) # most common(uncharacterized) / total found
                enr = (len(a)/float(len(item))) # unique / total found
                print (" **** Estimation:\nCoverage: ", ("{0:.2f}".format(cov)))#
                print ("Enrichment", ("{0:.2f}".format(enr)))
        print ("\n", a.most_common(25)) # 2088
        for j in sorted(dict(a.most_common(150)), key=lambda x :dict(a.most_common(150))[x]): # instead of 25 use all
            #print j, "\t", a[j]
            jfile.write(str(j)+ "\t"+ str( a[j])+"\n")
            
        jfile.write("\n")
        count += 1
        if count == len(fields):
            continue
    

    
    jfile.close()
    theInputFiles.append(jfile.name)
    print (jfile.name ," saved!")
    print ("************** done! ****************")

    if not with_fields:
       return protdata
    else:
        return protdata, fields


#readUniProtTab_sum("2pep_246_down.txt", n = 0, with_fields = False)
#"260_prot_merged_filtered_line.txt" ,"_bmq_ipr_go.txt","2pep_246_down.txt"


#teams_list = ["Man Utd", "Man City", "T Hotspur"]
#data = np.array([[1, 2, 1],
#                 [0, 1, 0],
 #                [2, 4, 2]])



#################################################################3##############
#### ENRICHMENT TEST    --------    PVALUE -   FDR ###########################
#########################################################3###############

# 9 May 2015
# http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python

def correct_pvalues_for_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):                
    """                                                                                                   
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """
    from numpy import array, empty                                                                        
    pvalues = array(pvalues) 
    n = float(pvalues.shape[0])                                                                           
    new_pvalues = empty(n)
    fdr_test = {}
    if correction_type == "Bonferroni":                                                                   
        new_pvalues = n * pvalues
    elif correction_type == "Bonferroni-Holm":                                                            
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        for rank, vals in enumerate(values):                                                              
            pvalue, i = vals
            new_pvalues[i] = (n-rank) * pvalue                                                            
    elif correction_type == "Benjamini-Hochberg":                                                         
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        values.reverse()                                                                                  
        new_values = []
        for i, vals in enumerate(values):                                                                 
            rank = n - i
            pvalue, index = vals                                                                          
            new_values.append((n/rank) * pvalue)                                                          
        for i in xrange(0, int(n)-1):  
            if new_values[i] < new_values[i+1]:                                                           
                new_values[i+1] = new_values[i]                                                           
        for i, vals in enumerate(values):
            #print vals
            pvalue, index = vals
            #print vals, new_pvalues
            new_pvalues[index] = new_values[i]                                                                                                                  
    return new_pvalues







def myformat(x):
    return ('%.11f' % x).rstrip('0').rstrip('.')

import numpy

def calc_benjamini_hochberg_corrections(p_values, num_total_tests):
    """
    Calculates the Benjamini-Hochberg correction for multiple hypothesis
    testing from a list of p-values *sorted in ascending order*.

    See
    http://en.wikipedia.org/wiki/False_discovery_rate#Independent_tests
    for more detail on the theory behind the correction.

    **NOTE:** This is a generator, not a function. It will yield values
    until all calculations have completed.

    :Parameters:
    - `p_values`: a list or iterable of p-values sorted in ascending
      order
    - `num_total_tests`: the total number of tests (p-values)

    """
    #cfile = open(cf, "r")
    
    bh_values = 0
    ph_values = 0
    prev_bh_value = 0
    print ("p value \t fdr")
    fdr_dict = {}
    
    for i, p_value in enumerate(p_values):
        #print i, p_value
        
        
        bh_value = float(p_value) * int(num_total_tests) / (int(i) + 1)
        #print float(bh_value)
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh_value = min(float(bh_value), 1)
        #print float(bh_value)
        # To preserve monotonicity in the values, we take the
        # maximum of the previous value or this one, so that we
        # don't yield a value less than the previous.
        bh_value = max(bh_value, prev_bh_value)
        prev_bh_value = bh_value
        #print p_value, float(bh_value)
        bh_values += bh_value
        ph_values += p_value
        pval = myformat(p_value)
        #print pval
        bhval = myformat(bh_value)
        #print bhval
        #yield bh_value
        fdr_dict[pval] = bhval
        
        
    print ("sum p_values = %s " %ph_values )
    print ("sum bh_values = %s" %bh_values)
    
        

    return fdr_dict


###############################################################################3

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
    print (header22)
    for i in res2[0:10]: print i
    bFile = open("output/_bmq_ipr_go.txt", "wt")
    bFile.write(str(header22))
    for aLine in res2:
        bFile.writelines(aLine+"\n")
    bFile.close()
    print  ("""Biomart Query | Hits : """, str(len(res2)))
    print (bFile, " saved!")
    
    
def FileCheck(kf):
    try:
      open(kf, "r")
      return 1
    except IOError:
      print ("Error: File does not appear to exist.\nLet 's create <bmq_ipr_go.txt>")
      get_bmq_ipr2go()
      return 0

#FileCheck("output/_bmq_ipr_go.txt")
#get_bmq_ipr2go()

#result = FileCheck("testfile")
#print result

###############################################################################33
##### PLOT 2AXES - BAR + LINE
############################################################################3

import sys, csv ,operator
def open_and_sort(af):

    data = csv.reader(open('File.csv'),delimiter=',')
    sortedlist = sorted(data, key=operator.itemgetter(0))    # 0 specifies according to first column we want to sort

    #now write the sorte result into new CSV file
    with open("NewFile.csv", "wb") as f:
        fileWriter = csv.writer(f, delimiter=',')
        for row in sortedlist:
            fileWriter.writerow(row)



##########################333
            
import matplotlib.pyplot as plt

import pandas as pd

def plot2yaxes(df, fields = True):

    # create bar chart.
    # df[['sales','net_pft']].unstack('STK_ID').plot(kind='bar', use_index=True)

    
    df = pd.read_csv(df ,delimiter = "\t", header=True) # delim_whitespace=True,header=None)

    print( df)
    #dfile = open(df, "r")
    

    #df[['counts','counts_all']].unstack('GO_ID').plot(kind='bar', use_index=True)

    # create line chart:
    #df[['sales_gr','net_pft_gr']].plot(kind='line', use_index=True)

    # df[['oddsratio','pvalue']].plot(kind='line', use_index=True

    fig = plt.figure()
    ax = df[['counts','IPRID']].unstack('totalcounts').plot(kind='bar', use_index=True)
                                    
    ax2 = ax.twinx()
    ax2.plot(df[['pvalue','oddsratio']].values, linestyle='-', marker='o', linewidth=2.0)

    show()
    
    #ax2.plot(ax.get_xticks(),df[['pvalue','oddsratio']].values, linestyle='-', marker='o', linewidth=2.0)

    # Change label limit

    #ax.set_ylim((-10, 80.))

#df = "output/A4I8L8_entry_ipr_name_gene_go_path_local_exi_ipr_enrichments.txt"
#plot2yaxes(df, fields = True)


###########################################################

#def plotlinescatter(df):

#########################################################
#############################################################3
from pylab import figure, show, legend, ylabel
import winsound
def plot_2ylines(df):

    # http://thomas-cokelaer.info/blog/2012/04/481/

    from collections import OrderedDict
    from collections import defaultdict
    time_begin = datetime.datetime.fromtimestamp(time.time())

    # Open data to plot

    dfile = open(df, "r")
    first = dfile.readline()
    first = first.strip()
    fields = first.split("\t")
    print (fields)
    
    gos = []
    counts = []
    pvals = []
    lines = []
    lines_id = {}
    count2pv_dict = defaultdict(list)
    for i in dfile.readlines()[1:]:
        i = i.strip()
        # print i
        #if "p-value" in fields: # IPR NAME
        #    pos = fields.index('p-value')
        #print pos
        #    ipr_name = line[pos] if len(line) > 2 else "null"
        #    ipr_names1.append(ipr_name)
        items = i.split("\t")
        ipr = items[0]
        lines.append(items)
        # l = [[0, 1, 'f'], [4, 2, 't'], [9, 4, 'afsd']]
        # l.sort(key=lambda x: x[2])
        if len(items) == 8:
            xlabel = items[0]
            #gos.append(items[0])
            pv =items[-2]
            #pvals.append(items[-2])
            cts = items[2]
            #counts.append(items[2])
            #count2pv_dict[cts].append(pv) #= pv
            count2pv_dict[pv] = cts# append(pv)   USING DICT
            joint = "\t".join(items)
            lines.append(items)
            lines_id[ipr]= joint
        
            

    #floats1 = " ".join(pvals)
    #floats = map(float, floats1.split())

    #print len(counts)
    #print len(pvals)
    
    #pvals = sorted(pvals)
    
    #counts = sorted(counts)
    
    
    d_sorted_by_value = OrderedDict(sorted(count2pv_dict.items(), key=lambda x: x[0]))
    for k, v in d_sorted_by_value.items():
        #print "%s: %s" % (k, v)
        pvals.append(k)
        counts.append(v)

    from operator import itemgetter
    pvals2 = []
    counts2 = []
    lf = open(str(dfile.name[:-4])+"_ordered.txt","w")
    order_lines = sorted(lines, key=itemgetter(-2)) #lines.sort(key=lambda x: x[-2]) :error
    for i in order_lines: # list inside list
        #print i[-2]
        val = i[2] if len(i) == 8 else "2000"
        
        #print val
        new_line2 = "\t".join(i)
        lf.write(new_line2+"\n")
        if int(val) < 1000:
            pvals2.append(i[-2])
            #print i[2]
            counts2.append(i[2])
        else:
            
            print (i[2],i[-2], val)

    lf.close()
    #print dfile.name[:-4]
    dfile.close()
    #global theInputFiles
    #theInputFiles.append(lf.name)

    title = fields[0] +" vs. " +fields[-2] + "\n" + dfile.name
    # create the general figure
    fig1 = figure()
    fig1.suptitle(str(title))
    
    

    # and the first axes using subplot populated with data
    ax1 = fig1.add_subplot(111)
    line1 = ax1.plot(pvals2, 'o-', label='p-value')
    #ax1.yaxis.tick_right()
    #ax1.yaxis.set_label_position("right")
    ylabel("P-value")
    #legend(loc="upper left")

    # now, the second axes that shares the x-axis with the ax1

    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
    line2 = ax2.plot(counts2, 'xr-', label='counts')
    
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ylabel("Positive Counts")
    #legend([line1], [line2], loc="upper left")
    #legend(loc="upper right")
    #ax1.legend()
    #ax2.legend()
    ax1.grid()
    ax1.set_xlabel(str(fields[0]))
    #ax1.legend(loc=0)


    # for the legend, remember that we used two different axes so, we need
    # to build the legend manually
    #legend([line1, line2], ["p-value", "counts"], loc="upper left")
    #fig1.legend(loc="upper right") #
    #legend([line1, line2], ['Item 1', 'Item 2'],loc="upper left")
    #legend()
    #host.xaxis_date()
    show()
    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    print ("(complete!)")
        
    time_end = datetime.datetime.fromtimestamp(time.time())
    print("Time elapsed: ", str(time_end - time_begin))

#df = "output/A4I8L8_entry_ipr_name_gene_go_path_local_exi_ipr_enrichments.txt"
#df1 = "output/A4I8Y2_entry_ipr_name_gene_go_path_local_exi_ipr_enrichments.txt"
#df2 = "output/A4I8Y2_entry_ipr_name_gene_go_path_local_exi_go_enrichments.txt"

#df3 = "output/E9AGK1_259_bs_uniprot_ipr_enrichments.txt"#"output/A4I1D9_248_bs_biomart_go_enrichments.txt"
#df5 = "output/E9AGK1_259_bs_uniprot_go_enrichments.txt"#"output/A4I1D9_248_bs_biomart_ipr_enrichments.txt"7
#df5 = "output/A4HRR5_259_bs_uniprot_ipr_enrichments.txt"#"output/A4I1D9_248_bs_biomart_ipr_enrichments.txt"7


#plot_2ylines(df3)
#plot_2ylines(df5)
    
# A4I8Y2_entry_ipr_name_gene_go_path_local_exi_go_enrichments
#df6 = "output/A4I8K7_259_bs_uniprot_ipr_enrichments.txt"
#df6 = "output/A4HUZ8_259_bs_uniprot"

#df6 = "output/E9AHP1_259_bs_uniprot_ipr_enrichments.txt"
#df7 = "output/E9AHP1_259_bs_uniprot_go_enrichments.txt"
##plot_2ylines(df6)
#plot_2ylines(df7)
####################################################################################3
#####  ENRICHMENT ANALYSIS - IPR AND GO -
#####################################################################################3

import scipy as sp
import scipy.stats
import numpy as np


def enrichment_analysis_domains(af, bf, kf): # n = 1000

    """Using the most common dict (COLLECTIONS) file from last functions,
    we are able to build a gene set default of the organism
    and compare it with input gene set"""


    # Gets tax id and retrieves uni ids list from uniprot web
    # then enrichment analysis with list of accessions

    #from os.path import basename
    # now you can call it directly with basename
    #print basename(af)
    #str(af[:-4])+
    print ("choice 1 = ", af)
    print ("choice 2 = ", bf)
    print ("dict ipr2go = " , kf)
    cfile = str(af[:-4])+"_go_enrichments.txt"
    cf = open(cfile, "w+")
    dfile = str(af[:-4])+"_ipr_enrichments.txt"
    df = open(dfile, "w+")
    
    iprs1 = []
    gos1 = []
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
                    gon_dict[j] = namego
                    


    for i in kfile.readlines():
        i = i.strip()
        entry = i.split("\t")
        print entry
        if len(entry) > 3 :
            namego = entry[-2]
            print namego
            for j in entry:
                if j[:3] == "GO:":
                    print j
                    gon_dict[j] = namego
                
    with open(af,"r") as f1:
        for i in f1.readlines():
            #print i
            i = i.strip()
            for j in i.split("\t"):
                #print j
                if j[0:3] == "IPR":
                    #print index(j)
                    #print j.split("; ")
                    for i in j.split("; "):
                        iprs1.append(i)
                        if i == "IPR027417": # in i:
                            print i
                elif j[0:3] == "GO:":
                    #print j.split("; ")
                    for i in j.split("; "):
                        gos1.append(i)

    iprs2 = []
    gos2 = []
    with open(bf,"r") as f2:
        for i in f2.readlines():
            #print i
            i = i.strip()
            for j in i.split("\t"):
                #print j
                if j[0:3] == "IPR":
                    #print j.split("; ")
                    for i in j.split("; "):
                        iprs2.append(i)
                elif j[0:3] == "GO:":
                    #print j.split("; ")
                    for i in j.split("; "):
                        gos2.append(i)
                        #else: # create dicts for id names
                        #default_data['item3'] = 3
                        #Easy as py.
                        #Another possible solution:
                        #default_data.update({'item3': 3})
                        #print i

    print ("TOTAL IPRs1 : ",len(iprs1),"\nTOTAL IPRs2 : " , len(iprs2))
    print ("TOTAL GOs1 : ",len(gos1),"\nTOTAL GOs2 : " , len(gos2))
                        
    one = set(iprs1) #uniq_iprs1)
    two = set(iprs2) # uniq_iprs2)

    # IPR IDs
    #**set(['a', 'c', 'b'])**
    print ("\nAll Domains ", af, bf)
    print ("1 Union 2 = ", "\t", len(one.union(two)))
    #a = most.common(one)
    
    print ("\nDomains in common: ", af, bf)
    print ("1 Intersection 2 = ", "\t", len(one.intersection(two)))
    #print one.intersection(two)

    #**set(['a', 'c', 'b', 'd'])**
    print ("\nDomains ", af, bf)
    print ("1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - one.intersection(two)))
    #**set(['d'])**

    #a = [1, 2, 3, 4, 5]
    #b = [9, 8, 7, 6, 5]
    #i for i, j in zip(a, b) if i == j]

    
    # GO TERMS
    three = set(gos1) #uniq_iprs1)
    four = set(gos2) # uniq_iprs2)

    #**set(['a', 'c', 'b'])**
    print ("\nAll GOs ", af, bf)
    print ("2 Union 1 = ", "\t", len(four.union(three)))
    #print four.union(three)
    
    print ("\nGOs in common: ", af, bf)
    print ("2 Intersection 1 = ", "\t", len(four.intersection(three)))
    #print one.intersection(two)

    
    #**set(['a', 'c', 'b', 'd'])**
    print ("\nGOs ", af, bf)
    print ("2 Union 1 - 2 Intersection 1 = ","\t", len(four.union(three)  - four.intersection(three)))

    ### IPR SUM ###
    c  = Counter(iprs1)
    print ("\nTOTAL iprs1 found : ", len(iprs1)) # 6597
    print ("UNIQUE iprs1 : ", len(c), "\n", c.most_common(10)) # 2088
    #for i in dict(c.most_common()):
     #   cf.write(str(i)+"\t"+str(c[i])+"\n")
        # print i, c[i]

    d  = Counter(iprs2)
    print "\nTOTAL IPRs2 found : ", len(iprs2) # 6597
    print "UNIQUE IPRs2 : ", len(d), "\n", d.most_common(10) # 3606
    #for i in dict(d.most_common()):
       # cf.write(str(i)+"\t"+str(d[i])+"\n")
        #print i, d[i]
    d_dict = dict(d.most_common())


    #e  = Counter(ippr_names)
    #print "\nTOTAL IPRs NAMES found : ", len(ippr_names) # 6597
    #print "UNIQUE IPR NAMES : ", len(e), "\n", e.most_common(10) # 3606
    #for i in dict(e.most_common(10)):
     #   print i, e[i]

    ### GO SUM ###
    f  = Counter(gos1)
    print "\nTOTAL GO IDS 1 found : ", len(gos1) # 6597
    print "UNIQUE GO IDs 1: ", len(f), "\n", f.most_common(10) # 1183
    #for i in dict(f.most_common()):
        #cf.write(str(i)+"\t"+str(f[i])+"\n")
        # print i, f[i]

    g  = Counter(gos2)
    print "\nTOTAL GO 2 found : ", len(gos2)# 6597
    print "UNIQUE GO 2: ", len(g), "\n", g.most_common(10) # 1183
    g_dict = dict(g.most_common())

    #for i in g_dict:
        #cf.write(str(i)+"\t"+str(g[i])+"\n")
        #print i, g_dict[i]
        #oddsratio, pvalue = stats.fisher_exact([[1100,6848], [11860,75292]])

    ### PVALUE + ODDSRATIO
    
    #print iprn_dict["IPR030470"]
    #iprn_dict["IPR030470"] = "null"
    df.write("IPR ID\tIPR Name\tcounts\ttotal counts\tquery\tall\tp-value\todds ratio\n") # lack id name. how to? fromgooripr_insert bmq_ipr2go
    results = []
    ipr_pv_dict = {}
    pvalues_ipr = []
    for j in sorted(dict(c.most_common()), key=lambda x :dict(c.most_common())[x]):
        #print j, f[j]
        if j in d_dict.keys():
            #print j, g_dict[j]
            n1 = c[j]
            n2 = d_dict[j]
            n3 = len(c)- n1
            n4 = len(d_dict) - n2
            
            oddsratio, pvalue = stats.fisher_exact([[n1,n2], [n3,n4]])
            #if n1 and n2 and n3 and n4 > 0 else "null"
            
            if j not in iprn_dict.keys():
                iprn_dict[j] = "null"
            pvalues_ipr.append(pvalue)
            ipr_pv_dict[j]=pvalue
            
            #pvalue = lambda pvalue: float(pvalue.replace('.',','))
            #oddsratio = oddsratio = lambda oddsratio: float(oddsratio.replace('.',',')) 
            
            #print j, n1 , n2, n3, n4,"\tp-value = ", pvalue, "\tOdds Ratio = ", oddsratio
            results.append(str(pvalue)+"\t"+str(oddsratio))
            df.write(str(j)+"\t"+str(iprn_dict[j])+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)+"\t"+str(pvalue)+"\t"+str(oddsratio)+"\n")

    cf.write("GO ID\tGO Term\tcounts\ttotal counts\tquery\tall\tp-value\todds ratio\n") # lack id name. how to? fromgooripr_insert bmq_ipr2go
    results = []
    go_pv_dict = {}
    pvalues_go = []
    for j in sorted(dict(f.most_common()), key=lambda x :dict(f.most_common())[x]):
        #print j, f[j]
        if j in g_dict.keys():
            #print j, g_dict[j]
            n1 = f[j]
            n2 = g_dict[j]
            n3 = len(f)- n1
            n4 = len(g_dict) - n2
            oddsratio, pvalue = stats.fisher_exact([[n1,n2], [n3,n4]])

            #result_arr = fisherExact(np.array([[5, 0], [1, 4]]))
            #print result_arr
            
            # SCIPY pvalue, uses hypergeometric distribution (not efficient)
            #x = [[n1,n2],[n3,n4]]
            #result_pv = sp.stats.fisher_exact(x)
            # print result_pv[1], "pv"
            #pvalue2 = result_pv[1] # for i in result_pv:
            #odds2 = result_pv[0]


            if j not in gon_dict.keys():
                gon_dict[j] = "null"
                #break

            #pval_go = myformat(pvalue)
            pvalues_go.append(pvalue)
            go_pv_dict[j] = pvalue

            #pvalue = lambda pvalue: float(pvalue.replace('.',',')) # convertfloat = lambda x: float(str(x).replace(',','.')) 
            #oddsratio = lambda oddsratio: float(oddsratio.replace('.',',')) # convertfloat = lambda x: float(str(x).replace(',','.')) 
            
            #print j, n1 , n2, n3, n4,"\tp-value = ", pvalue, "\tOdds Ratio = ", oddsratio
            results.append(str(pvalue)+"\t"+str(oddsratio))
            cf.write(str(j)+"\t"+str(gon_dict[j])+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)+"\t"+str(pvalue)+"\t"+str(oddsratio)+"\n") # "\t"+str(pvalue2)+"\t"+str(odds2)+
    
        
    
    #"oddsratio, pvalue = stats.fisher_exact([[1100,6848], [11860,75292]])
    #print pvalue, oddsratio
    #obs = array([[ 1100,  6848], [11860, 75292]])
    #from scipy.stats.contingency import expected_freq
    #expected_freq(obs)
    #array([[  1083.13438486,   6864.86561514],
    #   [ 11876.86561514,  75275.13438486]])

    #from scipy import stats
    #import numpy as np
    #obs = np.array(
    #    [[1100,6848],
    #   [11860,75292]])
    # stats.chi2_contingency(obs)

    ####################################       Benjamini Hochberg  : p-value adjustment                  #################
    
    pv_iprs = sorted(pvalues_ipr)
    fdr_ipr_dict = {}
    for i, p_value in enumerate(pv_iprs):
        fdr_ipr_dict[i] = p_value

    num_total_tests_ipr = len(pv_iprs)
    #benjamin_values_ipr = calc_benjamini_hochberg_corrections(pv_iprs, num_total_tests_ipr)

    bh_values_ipr = correct_pvalues_for_multiple_testing(pv_iprs, correction_type = "Benjamini-Hochberg")
    #print bh_values_ipr[1]
    pvalues_iprs = np.array(pv_iprs)

    

    new_ipr_array = np.dstack((pvalues_iprs,bh_values_ipr))
    #print new_array
    print "pval \t bh value"
    #for i,j in new_ipr_array[0]:
        #print i, j
   
    
     
    pv_gos = sorted(pvalues_go)
    fdr_go_dict = {}
    for i, p_value in enumerate(pv_gos):
        fdr_go_dict[i] = p_value
        
    num_total_tests_go = len(pv_gos)
    benjamin_values_go = calc_benjamini_hochberg_corrections(pv_gos, num_total_tests_go)
    #print benjamin_values_go.keys()

    bh_values_go = correct_pvalues_for_multiple_testing(pv_gos, correction_type = "Benjamini-Hochberg")
    #print bh_values_go[1]
    pvalues_gos = np.array(pv_gos)

    new_go_array = np.dstack((pvalues_gos,bh_values_go))
    #print new_array
    #print "pval \t bh value"
    # import numpy as np
    # data = np.genfromtxt("yourfile.dat",delimiter="\n")
    
    #with open(dfile) as ddf:
    ef = open(str(cfile[:-4])+"_adjusted.txt","w")
    #line1 = cf.readline()
    #print line1
    #line1 = line1.strip()
    #ef.write(str(line1)+"\n")
        
    for line in cf.readlines()[1:]:
        #lines = ddf.readlines()
        line = line.strip()
        
        #plot2yaxes(df, fields = True)

        items = line.split("\t")
        #pval = items[-2] if len(items)>3 else "null"
        #print pval
        for i,j in new_go_array[0]:
            #for key,value in go_pv_dict.iteritems():
            if i in items:
                print items, i, j
                ef.write("\t".join(items)+str(i+str(j))+"\n")
                #        go_pv_dict.pop(key,None) #del go_pv_dict[key]
                #continue
        
                #if i in go_pv_dict.values():
                # print i, j#, go_pv_dict[i].name
                # print "the key name is", mydic['key_name'].name_the_key(),
                # "and its value is", mydic['key_name']
                # for key, value in mydic.iteritems() :
                # print key, value
                # for line in file:
                # column.append(float(line[1:].split("\t")[3]))
                # column.sort()
        
        
    
    for line in cf.readlines():
        #print line
        line = line.strip()
        items = line.split("\t")
        p_val = items[6] if len(items) > 6 else "null"
        #print p_val
        for i in items:
            #print i
            
            if i in fdr_go_dict.values():
                print i #, benjamin_values_go[i]
                new_line = line + bh_values_go[1] # benjamin_values_go[i]
                print new_line
                #cf.write(mnew_line+"\n")
            elif i in bh_values_go:
                print i

    # Plot 2y Scatter plot : ID vs p-value
    
    plot_2ylines(dfile)
    plot_2ylines(cfile)
    
                
                
    

    

    # False Discovery Rate of the p-value 
        

    print cf, " saved!"
    print df, " saved!"
    print ef, " saved!"
    ef.close()
    cf.close()
    df.close()
    kfile.close()
    #name1 = os.path.basename(cf)
    #name2 = os.path.basename(df)

    #global theInputFiles
    #theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)
    Freq = 1000 # Set Frequency To 2500 Hertz
    Dur = 100 # Set Duration To 1000 ms == 1 second
    winsound.Beep(Freq,Dur)
    print "(complete!)"
    



#enrichment_analysis_domains("tests/2pep_246_down.txt", "tests/2peps.txt", n = 1000)
#enrichment_analysis_domains("tests/2pep_260.txt", "tests/2peps.txt", n = 1000)

#enrichment_analysis_domains("output/A4I2Z1_entry_ipr_name_gene_go_path_local_exi.txt", "output/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt", "output/_bmq_ipr_go.txt" ) #,n = 1000) # "tests/_bmq_ipr_go.txt"
#enrichment_analysis_domains("output/A4I2Z1_entry_ipr_name_gene_go_path_local_exi.txt", "output/A4I800_entry_ipr_name_gene_go_path_local_exi.txt", "output/_bmq_ipr_go.txt" ) #,n = 1000) # "tests/_bmq_ipr_go.txt"

#enrichment_analysis_domains("output/A4HRR5_259_bs_uniprot.txt", "output/A4HZ66_3990_bs_uniprot.txt", "output/_bmq_ipr_go.txt" )

# 260 vs all
# 246 vs all
#output/A4IBW8_entry_ipr_name_gene_go_path_local_exi.txt output/A4HU44_entry_ipr_name_gene_go_path_local_exi.txt

"""
output/A4IBW8_entry_ipr_name_gene_go_path_local_exi.txt  New Selection 1 !
choice 1 =  output/A4IBW8_entry_ipr_name_gene_go_path_local_exi.txt
choice 2 =  output/A4I2I6_entry_ipr_name_gene_go_path_local_exi.txt
dict ipr2go =  output/_bmq_ipr_go.txt
TOTAL IPRs1 :  610 
TOTAL IPRs2 :  11239
TOTAL GOs1 :  252 
TOTAL GOs2 :  4971"""
#################################################################3
##############################################################3
#######################################################3





def listCompare(list1,list2):

    #list1 = [1, 2, 3, 4, 5]
    #list2 = [5, 6, 7, 8, 9]
    if [item for item in list1 if item in list2]:
        print("Number was found")

    else:
        print("Number not in list")

#########################################################################
#############################    IPR VENN2   ############################
#########################################################################

def read_tab_enrichment_er_ipr(theFile1, theFile2):

    # Reads file as tab
    # Finds IPR
    # Finds True and False s
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
                        print i, " IPR FROM UNIPROT"
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
                        print i, " IPR FROM UNIPROT"
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
    
    #s = ['a','b','c']
    #f = ['a','b','d','c']
    one = set(uniq_iprs1)
    two =set(uniq_iprs2)
    
    print "\nDomains in common: ", af.name, bf.name
    print "1 Intersection 2 = ", "\t", len(one.intersection(two))
    #print one.intersection(two)

    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", af.name, bf.name
    print "1 Union 2 = ", "\t", len(one.union(two))
    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", af.name, bf.name
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - one.intersection(two))
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
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(two.union(one)  - two.intersection(one))
    
            

    
#read_tab_enrichment_ipr("output/5671_bmq_ac_ipr_go.txt", "output/5661_bmq_ac_ipr_go.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")
#read_tab_enrichment_ipr("outputs/5671_bmq_ac_ipr_go.txt", "outputs/5671_bmq_ac_ipr_go_short.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")
#read_tab_enrichment_ipr("outputs/5671_bmq_ac_ipr_go.txt", "outputs/5671_bmq_ac_ipr_go_short.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")

########################################################

def get_domain_venn_diagram(theFilename1, theFilename2, theFilename3, theFilename4, theOutput):
    "Entry  Entry name  Protein names   InterPro    Gene ontology IDs   Gene ontology (GO)  Gene names  Pathway Protein existence   Organism    Organism ID"
    "protein_accession	protein_name	entry_ac	entry_name	entry type	go_id	go_term	go_root	protein_database_name	tax_name	tax_id"
    # Given this format will produce countings for all InterPro entries and GO terms
    aFile = open(theFilename1, "rU")
    bFile = open(theFilename2, "rU")
    cFile = open(theFilename3, "rU")
    dFile = open(theFilename4, "rU") # Human
    eFile = open(theOutput, "wt")
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
        if len(aProtein_entry[3])== 0:
            continue
        elif len(aProtein_entry[3])== 9:
            Lin_domains[aProtein_entry[3]] = aRow
            if aProtein_entry[3] not in Lin_unique_domains:
                Lin_unique_domains.append(aProtein_entry[3])
        else:
            iprs1 = aProtein_entry[3]
            for i in iprs1.split("; "):
                #print i
                Lin_domains[i] = aRow
                if i not in Lin_unique_domains:
                    Lin_unique_domains.append(i)
            
    
    for aRow in bFile.readlines():
        #print aRow
        aProtein_entry = aRow.strip().split("\t")
        if len(aProtein_entry[3])== 0:
            continue
        elif len(aProtein_entry[3])== 9:
            Lma_domains[aProtein_entry[3]] = aRow
            if aProtein_entry[3] not in Lma_unique_domains:
                Lma_unique_domains.append(aProtein_entry[3])
        else:
            iprs2 = aProtein_entry[3]
            for i in iprs2.split("; "):
                #print i
                Lma_domains[i] = aRow
                if i not in Lma_unique_domains:
                    Lma_unique_domains.append(i)
    
    for aRow in cFile.readlines():
        #print aRow
        aProtein_entry = aRow.strip().split("\t")
        if len(aProtein_entry[3])== 0:
            continue
        elif len(aProtein_entry[3])== 9:
            Ldo_domains[aProtein_entry[3]] = aRow
            if aProtein_entry[3] not in Ldo_unique_domains:
                Ldo_unique_domains.append(aProtein_entry[3])
        else:
            iprs3 = aProtein_entry[3]
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
    print len(Lin_unique_domains), "\t in ", aFile.name  #1561
    print len(Lma_unique_domains), "\t in ", bFile.name #1570
    print len(Ldo_unique_domains), "\t in ", cFile.name  #1069
    
    for aRow in dFile.readlines():
        #print aRow
        aProtein_entry = aRow.strip().split("\t")
        if len(aProtein_entry[3])== 0:
            continue
        elif len(aProtein_entry[3])== 9:
            Hum_domains[aProtein_entry[3]] = aRow
            if aProtein_entry[3] not in Hum_unique_domains:
                Hum_unique_domains.append(aProtein_entry[3])
        else:
            iprs3 = aProtein_entry[3]
            for i in iprs3.split("; "):
                #print i
                Hum_domains[i] = aRow
                if i not in Hum_unique_domains:
                    Hum_unique_domains.append(i)
    print len(Hum_unique_domains), "\t in ", dFile.name  #1069



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
    four = set(Hum_unique_domains)
    
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



    
    ####### Human ####
    print "\nDomains in common: ", aFile.name, dFile.name
    print "2 Intersection 1 = ", "\t", len(one.intersection(four))
    
    look = list(one.intersection(two))
    print look[0:25]

    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", aFile.name, dFile.name
    print "2 Union 1 = ", "\t", len(one.union(four))
    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", aFile.name, dFile.name
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(one.union(four)  - one.intersection(four))

    



       

#get_domain_venn_diagram ("output/venns/5664_ac_entry_ipr_name_gene_go_path_local_exi.txt",
 #                        "output/venns/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt",
  #                      "output/venns/5661_ac_entry_ipr_name_gene_go_path_local_exi.txt",
   #                      "output/venns/9606_ac_entry_ipr_name_gene_go_path_local_exi.txt", 
    #                    "output/venns/Lin_complex_venn.txt")

#get_domain_venn_diagram ("output/venns/5661_ac_entry_ipr_name_gene_go_path_local_exi.txt",
 #                        "output/venns/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt",
  #                      "output/venns/5664_ac_entry_ipr_name_gene_go_path_local_exi.txt",
   #                      "output/venns/9606_ac_entry_ipr_name_gene_go_path_local_exi.txt", 
    #                    "output/venns/Lin_complex_venn.txt")

#get_domain_venn_diagram ("output/venns/5664_ac_entry_ipr_name_gene_go_path_local_exi.txt",
 #                        "output/venns/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt",
  #                      "output/venns/5661_ac_entry_ipr_name_gene_go_path_local_exi.txt",
   #                      "output/venns/9606_ac_entry_ipr_name_gene_go_path_local_exi.txt", 
    #                    "output/venns/Lin_complex_venn.txt")

#get_domain_venn_diagram ("output/venns/5664_ac_entry_ipr_name_gene_go_path_local_exi.txt",
 #                        "output/venns/5661_ac_entry_ipr_name_gene_go_path_local_exi.txt",
  #                      "output/venns/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt",
   #                      "output/venns/9606_ac_entry_ipr_name_gene_go_path_local_exi.txt", 
    #                    "output/venns/Lin_complex_venn.txt")

#a = [1, 2, 3, 4, 5]
    #b = [9, 8, 7, 6, 5]
    #i for i, j in zip(a, b) if i == j]
##########################################################
    
theFilenames = ["_bmq_ipr_go.txt", "outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt", "outputs/5671_bmq_ac_ipr_go.txt"]
                #"_bmq_ipr_go.txt",
                #"",
                #"",
                #""]


#print 'Reading UniProt tab file -----------------------------------------------'
#filename = 'uniprot-organism_L_infantum.tab'#"output/5671_bmq_ac_ipr_go.txt" # 
#prot_data, fields = readUniProtTab(filename, with_fields=True)
#print 'there are', len(prot_data), 'proteins in file', filename
#print

# import summary_functions_tk
#for i in theFilenames:
 #   print 'Reading UniProt tab file -----------------------------------------------'
  #  filename = i
   # prot_data, fields = readUniProtTab1(filename, with_fields=True)
    #print 'there are', len(prot_data), 'proteins in file', filename



#####

def read_Ids_tab_delimited (filename):
    ids = []
    a = open(filename)
    for linha in a:
        linha = linha.strip()
        id = linha.split('\t',1)[0]
        ids.append(id.strip())
    a.close()
    return ids

#print 'Reading ids ----------------------------------------------------'
#ids_filename1 = 'ColIds_2peps.txt'
#ids_filename2 = 'ColIds_all.txt'
#ids1 = read_Ids_tab_delimited(ids_filename1)
#ids2 = read_Ids_tab_delimited(ids_filename2)

#############3



def load_proteins_accessions(theFilename):
    aProteins = {}
    aFile = open(theFilename, "rt")
    for aLine in aFile.readlines()[1:]:
        aRow = aLine.strip().split("\t")# FORMAT = 5671_ac_entry_ipr_name_gene_go_goid_path_local_exi_tax
        print aRow
        if len(aRow)== 11:
            aProteins[str(aRow[0])].append({'name':aRow[1],'uniprot':aRow[0], 'interpro':aRow[3], 'gene_name':aRow[5], 'go_term':aRow[6],'protein_name':aRow[4]})
        else:
            aProteins[str(aRow[0])] = [{'name':aRow[1],'uniprot':aRow[0],'interpro':aRow[3], 'gene_name':aRow[5], 'go_term':aRow[6],'protein_name':aRow[4]}]
    aFile.close()
    return aProteins




def count_ipr_from_table_list(theFilename):
    aFile = open(theFilename, "rU")
    #bFile = open (theFilename2, "wt")
    #cFile = open(theFilename3, "wt")
    UniProt_accessions = []
    uniprot_ac = []
    interpro_ac = []
    count_ipr = 0
    ac_ipr = 0
    InterPro_accessions = []
    Protein_entries = []
    aProtein_ac_ipr =[]
    for aProtein in aFile:
        #print aProtein
        #aProtein = aProtein.split("\t")
        org_tax = aProtein.split("\t")[-1]
        org_name = aProtein.split("\t")[-2]
        uniprot_ac.append(aProtein[2])
        aProtein_ac = aProtein.split("\t")[0]
        aProtein_ipr = aProtein.split("\t")[3]
        aProtein_ac_ipr.append(aProtein_ac)
        aProtein_ac_ipr.append(aProtein_ipr) # ['A4HY57', 'IPR013785; IPR006638; IPR003698; IPR007197; ']
        #print aProtein_ac_ipr
        aProtein_name = aProtein.split("\t")[2]
        for i in aProtein_ac_ipr:
            #print i
            if i[0:3] == "IPR":
                ac_ipr = ac_ipr + 1
                ipr = i.split("; ") #['IPR001816', 'IPR014039', 'IPR009060', '']
                for i in ipr:
                    interpro_ac.append(i)
                    count_ipr = count_ipr + 1
        del aProtein_ac_ipr[:]
    print "\n\nOrganism :", org_name
    print "Proteins Found in UniProt BioServices:" , len(uniprot_ac)
    print "Total InterPro accessions :",len(interpro_ac), count_ipr #
    InterPro_acs = []
    for ipr in interpro_ac:
        if ipr not in InterPro_acs:
            InterPro_acs.append(ipr)
            #cFile.writelines(aProtein+"\n")
    print "Unique Interpro accessions :", len(InterPro_acs) #
    print "From total UniProt proteins ", len(uniprot_ac), "we found IPR annotations for ",ac_ipr," proteins."
    print InterPro_acs[0:10]
    #bFile.close()
    #cFile.close()
    return uniprot_ac
    

def count_go_terms_from_table_list(theFilename):
    aFile = open(theFilename, "rU")
    #bFile = open(TheOutput1, "wt")
    #cFile = open(TheOutput2, "wt")
    interpro_ac =[]
    uniprot_ac = []
    goterms_ac = []
    iprs_ac = []
    go_term_acs = []
    count_go_term = 0
    for aLine in aFile:
        aLine = aLine.split("\t")
        org_name = aLine[-2]
        uniprot_ac.append(aLine[2])
        for item in aLine:
            if item[0:3] == "IPR":
                iprs = item.split("; ")
                for i in iprs:
                    iprs_ac.append(i)
                    interpro_ac.append(aLine[3])
        for item in aLine:
            if item[0:3] == "GO:":
                #print item
                go_terms = item.split("; ") # ['GO:0005840', 'GO:0003735', 'GO:0006412']
                for i in go_terms:
                    go_term_acs.append(i)
                count_go_term = count_go_term + 1
    print "\n\nOrganism :", org_name
    print "Proteins Found in UniProt BioServices:" , len(uniprot_ac)
    print "InterPro accessions : %i" % (len(interpro_ac)), len(iprs_ac)# 53100
    print "GO Terms : %s " % (len(go_term_acs)) # 8892
    print "Proteins with go terms:", count_go_term
    unique_go = []
    for go in go_term_acs:
        if go not in unique_go:
            unique_go.append(go)
    print "Unique GO terms: ", len(unique_go)
    
def count_path_from_table_list(theFilename):
    aFile = open(theFilename, "rU")
    #bFile = open (theFilename2, "wt")
    #cFile = open(theFilename3, "wt")
    count_no_path = 0
    count_path = 0
    path = []
    for aProtein in aFile:
        #print aProtein
        org_name = aProtein.split("\t")[-2]
        aProtein_ac = aProtein.split("\t")[0]
        aProtein_ipr = aProtein.split("\t")[3]
        aProtein_ac_ipr = aProtein_ac+ aProtein_ipr
        aProtein_names = aProtein.split("\t")[2]
        aProtein_path = aProtein.split("\t")[7]
        #print aProtein_path
        if aProtein_path =="Pathway":
            print "Pathway"
        elif aProtein_path == "":
            #print "Pathway not Found in ", aProtein_ac
            count_no_path = count_no_path + 1
        else:
            #print "\n",aProtein_ac
            #print aProtein_names
            #print aProtein_path
            count_path = count_path + 1
            path.append(aProtein_path.split("; ")[0])
            #for i in path:
                #print aProtein_ac, i
    for i in Counter(path):
        print i
    print "Pathways found:"
    print "Proteins found with pathways: ",count_path
    print "Proteins found with no pathways: ",count_no_path
    print count_no_path + count_path
    print Counter(path)
    #aProtein_path = aProtein_path.split("; ")
    #print aProtein_path
            #print i
                    #count_ipr = count_ipr + 1
    #print "Proteins Found in UniProt Bioservices:" , len(uniprot_ac)


def count_gene_from_table_list(theFilename):
    aFile = open(theFilename, "rU")
    #bFile = open (theFilename2, "wt")
    #cFile = open(theFilename3, "wt")
    count_gene = 0
    gene = []
    aGene_Dict = {}
    for aProtein in aFile:
        #print aProtein
        org_name = aProtein.split("\t")[-2]
        aProtein_ac = aProtein.split("\t")[0]
        aProtein_ipr = aProtein.split("\t")[3]
        aProtein_ac_ipr = aProtein_ac + aProtein_ipr
        aProtein_names = aProtein.split("\t")[2]
        aProtein_gene = aProtein.split("\t")[6]
        aGene_Dict[aProtein_ac]= aProtein_gene
    #for i in aGene_Dict:
        #print i
    #for i in aProtein_gene.split(" "):
    #print aProtein_ac, i
    print "Genes Found: ",len(aGene_Dict)
            

def count_protein_existence(TheFilename):
    aFile = open(TheFilename, "rU")
    #bFile = open("output/apendix.txt", "wt")
    aLines = []
    status = []
    count_plevel = 0
    count_tlevel= 0
    count_inferred = 0
    count_predicted = 0
    count_uncertain = 0
    """1. Evidence at protein level
    2. Evidence at transcript level
    3. Inferred from homology
    4. Predicted
    5. Uncertain"""
    for aLine in aFile.readlines():
        aProtein = aLine.split("\t")
        org_tax = aLine.split("\t")[-1]
        org_name = aLine.split("\t")[-2]
        aProtein_ac = aLine.split("\t")[0]
        aProtein_ipr = aLine.split("\t")[3]
        #aProtein_genes = aLine.split("\t")[6]
        aProtein_status = aLine.split("\t")[8]
        #print aProtein_genes
        #print aProtein_ac, aProtein_ipr
        if aProtein_status not in status:
            status.append(aProtein_status)
        if  "Evidence at protein level" in aProtein_status:
            count_plevel = count_plevel + 1
        elif "Evidence at transcript level" in aProtein_status:
            count_tlevel = count_tlevel + 1
        elif "Inferred from homology" in aProtein_status:
            count_inferred = count_inferred + 1
        elif "Predicted" in aProtein_status:
            count_predicted = count_predicted +1
        elif "Uncertain" in aProtein_status:
            count_uncertain = count_uncertain + 1
        else:
            print  "Found other status:", aProtein_status
        aLines.append(aProtein_status)
    print "BioServices - UniProt KB\n"
    print org_name, org_tax
    print "Proteins found:", len(aLines)
    #print aLines[0:8]
    print status[0]
    print count_plevel, "Evidence at protein level"#9
    print count_tlevel, "Evidence at transcript level"#36
    print count_inferred, "Inferred from homology"#933
    print count_predicted, "Predicted"#7316
    print count_uncertain, "Uncertain"#0>Total=8294
    aFile.close()
    print "done!"


#################################################################################3    
##############################################################################
### INTERPRO - BIOMART: INPUT FILE <ac_ipr>
#############################################################################
#INPUT: File_Format<

def load_GO_list(theFilename):
    """
    !date: 2014/05/09 09:05:20
    !Mapping of InterPro entries to GO
    !
    InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:DNA binding ; GO:0003677"""
    aFile = open(theFilename, "rU")
    aDomainDict = {}
    aLines = []
    for aRow in aFile.readlines():
        aDomain = aRow.strip()
        #print aDomain[9:19]
        if aDomain[0] == "!":
            print aDomain[1:]
        else:
            aLines.append(aDomain)
            aDomainDict[aDomain[9:19]] = aDomain[19:]
    aFile.close()    
    print "There are ", len(aLines), " domains in InterPro DB with GO Annotations." # 27857
    aFile.close()
    return aDomainDict

def get_interpro2go(theFileName, aDomainDict):
    """
    ['A0AAK3', 'IPR023186']['A0AAK3', 'IPR001910']['A1Y2C6', 'IPR002085']
    """
    aFile = open(theFileName, "rU")
    #bFile = open(theOutput, "wt")
    accessions = []
    domains =[]
    Lines = []
    ac2ipr2go = {}    
    for aRow in aFile.readlines():
        aProtein_entry = aRow.strip().split("['")
        #print aProtein_entry
        Lines.append(aRow)
        for i in aProtein_entry:
            if i[0:6] not in accessions: # See UniProt Pattern
                accessions.append(i[0:6])
                #print i[0:6]
            if i[10:19] not in domains:
                domains.append(i[10:19])
                #print i[10:19]
                #Lines.append(i)
            for key in aDomainDict:
                #print key
                if key == i[10:19]:
                    print i[0:6]
                    print i[10:19]
                    print key, aDomainDict[key]
    print " Lines", len(Lines)#14058
    print " Unique Accessions", len(accessions)#5350
    print " Unique Domains", len(domains)#3494
    #bFile.write(aLine)
    aFile.close()
    #bFile.close()



        
    
### Search for LEIIN orthologs
### Find Transporter
    ## find_term("catalase")\"IPR011614"\ "AC"
    ## uniprot_db_match()
    ##

##################################################################
def get_orthologues(theFileName, Lma_Proteins):
    aFile = open(theFileName, "rU")
    for i in aFile.readlines():
        print i
        if Lma_Proteins["gene_name"] in i:
            print i




##################################################################
#biopython
"""['A4I1F7',
('GO', 'GO:0019799', 'F:tubulin N-acetyltransferase activity', 'IEA:UniProtKB-HAMAP'),
('GO', 'GO:0071929', 'P:alpha-tubulin acetylation', 'IEA:InterPro'),
('GO', 'GO:0070507', 'P:regulation of microtubule cytoskeleton organization', 'IEA:UniProtKB-HAMAP'),
('InterPro', 'IPR016181', 'Acyl_CoA_acyltransferase'),
('InterPro', 'IPR007965', 'Alpha-TAT')]"""

"""5671 ASNA_LEIIN  A4HUY0 IPR025723; IPR016300; IPR027542; IPR027417;
ATPase ASNA1 homolog (EC 3.6.-.-) (Arsenical pump-driving ATPase homolog) (Arsenite-stimulated ATPase)
LinJ11.0710 LinJ_11_0710
GO:0005524; GO:0016887; GO:0005783; GO:0046872; GO:0045048; GO:0006810
ATP binding; ATPase activity; endoplasmic reticulum; metal ion binding; protein insertion into ER membrane; transport
Cytoplasm. Endoplasmic reticulum.   Inferred from homology  Leishmania infantum"""



def save_uniprot_ac(theProteins, theInput, theOutput):
    aFile = open(theInput, "rt")
    bFile = open(theOutput, "wt")
    acs =[]
    for aLine in aFile:
        aProtein_entry = aLine.strip().split('\t')
        for i in aProtein_entry:
            if len(i)== 6:
                acs.append(i)
    print len (acs)
    bFile.write(aLine)
    aFile.close()
    bFile.close()


################################################################################################3

theFilenames = ["output/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt",
                "output/5671_bmq_ac_ipr_go.txt",
                "",
                "",
                "",
                ""]


def from_bs_interpro_get_stats():

    """
    ac ipr go go_term
    """

    aFile = open("output/quickGO_interpro_db.goa", "rU")
    #bFile = open("output/biomartview_Interpro_ac.txt", "wt")
    #cFile = open("output/biomartview_GOterms_ac.txt", "wt")

    # ['UniProtKB', 'A1Y2C6', '', '', 'GO:0008270', 'GO_REF:0000002', 'IEA', 'InterPro:IPR002085|InterPro:IPR002364|InterPro:IPR013149', 'F', 'Zeta-crystallin/NADPH-oxidoreductase-like protein', 'A1Y2C6_LEIIN', 'protein', 'taxon:5671', '20140111', 'InterPro']

    count = 0
    
    interpro_ac =[]
    interpro_ac2 =[]
    uniprot_ac = []
    goterms_ac = []
    ac_ids = []
    tot_ac = []
    tot_go = []

    for aLine in aFile.readlines()[6:]:
        aLine = aLine.strip()
        #print aLine
        #aLine = aLine.split("\t")
        ac = aLine.split("\t")[1]
        tot_ac.append(ac)
        
        go_term = aLine.split("\t")[4]
        tot_go.append(go_term)
        
        iprs = aLine.split("\t")[7]
        go_root = aLine.split("\t")[8]
        go_name = aLine.split("\t")[9]

        if ac not in ac_ids:
            ac_ids.append(ac)

        iprss = iprs.split("|")
        for i in iprss:
            #print i
            if len(iprss)== 1:
                count = count +1
                interpro_ac.append(i)
            else:
                count = count + len(iprss)
                interpro_ac.append(i)
                #cFile.write(str(aLine[0:2])+str(aLine[6:9]))


    for i in interpro_ac:
        if i not in interpro_ac2:
            interpro_ac2.append(i)
            
    unique_goterms_ac = []
    for i in tot_go:
        if i not in unique_goterms_ac:
            unique_goterms_ac.append(i)


    print "Total ACS : ", len(tot_ac)          # 10716
    print "TOTAL ACS with  IPR: ", len(ac_ids) # 3488
    print "TOTAL IPR : ", len(interpro_ac)     # 16126
    print "Unique IPRS : ", len(interpro_ac2)  # 2123
    print "TOTAL GOs : ", len(tot_go)
    print "UNIQUE GOs : ", len(unique_goterms_ac)

    for i in interpro_ac[0:20]:
        print i

    print count #85136704, 33886

    # TASK - Filtering unique_ac

    

    print "Unique GO terms : ", len(unique_goterms_ac) # 1139, 0

    unique_interpro_ac = []
    for i in interpro_ac2:
        if i not in unique_interpro_ac:
            unique_interpro_ac.append(i)
            #bFile.write(str(i))

    print "Unique IPR : ", len(unique_interpro_ac) # 14058, 2123

    unique_uniprot_ac = []
    for i in uniprot_ac:
        if i not in unique_uniprot_ac:
            unique_uniprot_ac.append(i[0])
            print "Unique AC : ",len(ac_ids) # unique_uniprot_ac) # 5349

    print "DONE"



####
def summary_funtions_1(theFilenames):

    print "\nUniProt by Tax ID "
    print '[1] Loading protein accessions '
    print theFilenames[0]
    #Lin_Proteins = load_proteins_accessions(str(theFilenames[0])) # "output/5671_ac_entry_ipr_name_gene_go_path_local_exi_short.txt")

    print '[2] Summary Analysis '
    
    count_ipr_from_table_list(str(theFilenames[0]))
    count_go_terms_from_table_list(str(theFilenames[0]))
    count_path_from_table_list(str(theFilenames[0]))
    count_gene_from_table_list(str(theFilenames[0]))
    count_protein_existence(str(theFilenames[0]))

    print "\nInterPro Biomart by tax id "
    print '[3] Map InterPro 2 GO with BiomartQuery  '
    #get_interpro2go("output/5671_ac_ipr_interpro-bmq.txt", "output/interpro2go.txt")
    theDomains = load_GO_list("output/interpro2go.txt")
    get_interpro2go(str(theFilenames[1]), theDomains) # output/5671_ac_ipr_interpro-bmq.txt"
    
    ####################################################################
    
    print '\n[4] Map InterPro 2 GO with BiomartView '
    

#summary_funtions_1(theFilenames)








#####
##### BIOMARTVIEW SUMMARY   ####
#####

##### QUICK GO              ####


# TASK - read and count

def read_count_quickgo_bmv(theFilenameX):

    """
    ['UniProtKB', 'Protein', 'Accession',
    'InterPro', 'Entry', 'Accession',
    'InterPro', 'Entry', 'Name',
    'Source','Protein', 'Database',
    'NCBI', 'Taxonomy', 'ID',
    'Signature', 'Accession',
    'GO', 'ID',
    'GO', 'Term', 'Name',
    'GO', 'Root', 'Term', '(Process', '/', 'Component', '/', 'Function)']
    """

    aFile = open("output/quickGO_interpro_db.goa", "rU")
    #bFile = open("output/biomartview_Interpro_ac.txt", "wt")
    #cFile = open("output/biomartview_GOterms_ac.txt", "wt")

    # ['UniProtKB', 'A1Y2C6', '', '', 'GO:0008270', 'GO_REF:0000002', 'IEA', 'InterPro:IPR002085|InterPro:IPR002364|InterPro:IPR013149', 'F', 'Zeta-crystallin/NADPH-oxidoreductase-like protein', 'A1Y2C6_LEIIN', 'protein', 'taxon:5671', '20140111', 'InterPro']

    count = 0
    
    interpro_ac =[]
    interpro_ac2 =[]
    uniprot_ac = []
    goterms_ac = []
    ac_ids = []
    tot_ac = []
    tot_go = []

    for aLine in aFile.readlines()[6:]:
        aLine = aLine.strip()
        #print aLine
        #aLine = aLine.split("\t")
        ac = aLine.split("\t")[1]
        tot_ac.append(ac)
        
        go_term = aLine.split("\t")[4]
        tot_go.append(go_term)
        
        iprs = aLine.split("\t")[7]
        go_root = aLine.split("\t")[8]
        go_name = aLine.split("\t")[9]

        if ac not in ac_ids:
            ac_ids.append(ac)

        iprss = iprs.split("|")
        for i in iprss:
            #print i
            if len(iprss)== 1:
                count = count +1
                interpro_ac.append(i)
            else:
                count = count + len(iprss)
                interpro_ac.append(i)
                #cFile.write(str(aLine[0:2])+str(aLine[6:9]))


    for i in interpro_ac:
        if i not in interpro_ac2:
            interpro_ac2.append(i)
            
    unique_goterms_ac = []
    for i in tot_go:
        if i not in unique_goterms_ac:
            unique_goterms_ac.append(i)


    print "Total ACS : ", len(tot_ac)          # 10716
    print "TOTAL ACS with  IPR: ", len(ac_ids) # 3488
    print "TOTAL IPR : ", len(interpro_ac)     # 16126
    print "Unique IPRS : ", len(interpro_ac2)  # 2123
    print "TOTAL GOs : ", len(tot_go)
    print "UNIQUE GOs : ", len(unique_goterms_ac)

    for i in interpro_ac[0:20]:
        print i

    print count #85136704, 33886

    # TASK - Filtering unique_ac

    

    print "Unique GO terms : ", len(unique_goterms_ac) # 1139, 0

    unique_interpro_ac = []
    for i in interpro_ac2:
        if i not in unique_interpro_ac:
            unique_interpro_ac.append(i)
            #bFile.write(str(i))

    print "Unique IPR : ", len(unique_interpro_ac) # 14058, 2123

    unique_uniprot_ac = []
    for i in uniprot_ac:
        if i not in unique_uniprot_ac:
            unique_uniprot_ac.append(i[0])
            print "Unique AC : ",len(ac_ids) # unique_uniprot_ac) # 5349

    print "DONE"


# read_count_quickgo_bmv("output/quickGO_interpro_db.goa")
########################33



#for i in uni_ids:
def file2list(theFilename):
 
    
    print "File Formart :\nA4HT31','IPR008974','TRAF-like','UniProt/TrEMBL','Leishmania infantum','5671','IPR008974','GO:0005515','protein binding','function']"
    
    #aProteins = {}
    aFile = open(theFilename, "rt") # 'a', 'r+'
    #for aLine in aFile.readlines()[1:]:
     #   aRow = aLine.strip().split("\t")# FORMAT1 = ac_ipr_name_db_taxn_taxid_ipr_go_goid_root
      #  #aRow = aLine.strip().split("\t")# FORMAT2 = 5671_ac_entry_ipr_name_gene_go_goid_path_local_exi_tax
        #print aRow
       # if len(aRow)== 10: # list index out of range, key error, 
        #    aProteins[str(aRow[0])].append({'ipr_name':aRow[2],'uniprot':aRow[0], 'interpro':aRow[1], 'tax_id':aRow[5], 'go_id':aRow[7],'go_term':aRow[8], 'go_root':aRow[9]})
        #else:
         #   aProteins[str(aRow[0])] = [{'ipr_name':aRow[2],'uniprot':aRow[0],'interpro':aRow[1], 'tax_id':aRow[5], 'go_id':aRow[7],'go_term':aRow[8], 'go_root':aRow[9]}]
    #aFile.close()
    #return aProteins
    
    #af = open(theFilename, "rt")

    acs = []
    accs = []
    iprs = []
    ipprs = []
    ipr_names = []
    ippr_names = []
    go_ids = []
    goo_ids = []
    go_terms = []
    goo_terms = []
    go_roots = []
    for i in aFile.readlines()[1:]: # LIST OF LINES
        i = i.strip().split("\t")
        #print i
        ac = i[0]
        accs.append(ac)
        if ac not in acs:
            acs.append(ac)
            
        ipr = i[1]
        ipprs.append(ipr)
        if ipr not in iprs:
            iprs.append(ipr)
            
        ipr_name = i[2]
        ippr_names.append(ipr_name)
        if ipr_name not in ipr_names:
            ipr_names.append(ipr_name)
            
        go_id = i[7] if len(i) == 10 else 'null' # gotdata = dlist[1] if len(dlist) > 1 else 'null'
        goo_ids.append(go_id)
        if go_id not in go_ids:
            go_ids.append(go_id)

        go_term = i[8] if len(i) == 10 else 'null'
        goo_terms.append(go_term)
        if go_term not in go_terms:
            go_terms.append(go_term)

        go_root = i[9] if len(i) == 10 else 'null' 
        go_roots.append(go_root)

    c  = Counter(accs)
    print "\nTOTAL ACS found : ", len(accs) # 6597
    print "UNIQUE ACS : ", len(c), "\n", c.most_common(10) # 2088
    for i in dict(c.most_common(10)):
        print i, c[i]

    d  = Counter(ipprs)
    print "\nTOTAL IPRs found : ", len(ipprs) # 6597
    print "UNIQUE IPRs : ", len(d), "\n", d.most_common(10) # 3606
    for i in dict(d.most_common(10)):
        print i, d[i]

    e  = Counter(ippr_names)
    print "\nTOTAL IPRs NAMES found : ", len(ippr_names) # 6597
    print "UNIQUE IPR NAMES : ", len(e), "\n", e.most_common(10) # 3606
    for i in dict(e.most_common(10)):
        print i, e[i]
        
    f  = Counter(goo_ids)
    print "\nTOTAL GO IDS found : ", len(goo_ids) # 6597
    print "UNIQUE GO IDs : ", len(f), "\n", f.most_common(10) # 1183
    for i in dict(f.most_common(10)):
        print i, f[i]

    g  = Counter(goo_terms)
    print "\nTOTAL GO TERMS found : ", len(goo_terms)# 6597
    print "UNIQUE GO Terms : ", len(g), "\n", g.most_common(10) # 1183
    for i in dict(g.most_common(10)):
        print i, g[i]
        
    h  = Counter(go_roots)
    print "\nTOTAL GO ROOTS found : ", len(go_roots)
    print "UNIQUE GO ROOTS : ", len(h), "\n", h.most_common(10)
    
    # lines_of_text = ["a line of text", "another line of text", "a third line"]
    # fh.writelines(lines_of_text)
    # print 'The story of {0}, {1}, and {other}.'.format('Bill', 'Manfred', other='Georg')
    aFile.close()
        
        

"""#https://docs.python.org/2/library/collections.html
c = Counter(a=3, b=1)
>>> d = Counter(a=1, b=2)
>>> c + d                       # add two counters together:  c[x] + d[x]
Counter({'a': 4, 'b': 3})
>>> c - d                       # subtract (keeping only positive counts)
Counter({'a': 2})
>>> c & d                       # intersection:  min(c[x], d[x])
Counter({'a': 1, 'b': 1})
>>> c | d                       # union:  max(c[x], d[x])
Counter({'a': 3, 'b': 2})

d.rotate(-1)

https://docs.python.org/2/library/collections.html"""
###############uni_ids = file2list("output/5671_bmq_ac_ipr_go.txt")
#from_uni_id_get_bs_ac2ipr2go(uni_ids)
    
#########################################################################






def get_enrichment_with_bg(thebgfile, theinfile):

    # Reads file as tab
    # Finds IPR
    # Finds True and False s
    
    from collections import defaultdict
    
    af = open(thebgfile, "r")
    bf = open(theinfile, "r")

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

    uniq_gos1_ac = [] 
    
    ipr_dict1 = defaultdict(list)
    ipr_dict2 = defaultdict(list)
    go_dict1 = defaultdict(list)
    go_dict2 = defaultdict(list)
    
    for i in af.readlines():
        #print i
        
        i = i.strip()
        line = i.split("\t")
        ac = line[0]
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
                ipr_dict1[i] = line
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
                        print i, " IPR FROM UNIPROT"
                        iprs1.append(i)
                        if i not in uniq_iprs1:
                            uniq_iprs1.append(i)
                    
        for i in line:
            if i.startswith("GO:"): # GO
                go_dict1[i] = line
                go = i
                #print go
                gos1.append(i)
                if i not in uniq_gos1:
                    uniq_gos1_ac.append(ac)
                    uniq_gos1.append(i)

    uniq_gos2_ac = []
    for i in bf.readlines():
        #print i
        
        i = i.strip()
        line = i.split("\t")
        ac = line[0]
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
                ipr_dict2[i] = line
                ipr2 = i
                #print i
                if len (ipr2) == 9:
                    #print ipr1, # does: IPR2GO Biomart
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
                go_dict2[i] = line
                go = i
                #print go #
                gos2.append(go)
                if i not in uniq_gos2:
                    uniq_gos2_ac.append(ac)
                    
                    uniq_gos2.append(i)
    
    #s = ['a','b','c']
    #f = ['a','b','d','c']
    one = set(uniq_iprs1)
    two =set(uniq_iprs2)

    three = set(uniq_gos1)
    four = set(uniq_gos2)
    
    print "\nDomains in common: ", af.name, bf.name
    print "1 Intersection 2 = ", "\t", len(one.intersection(two))
    #print one.intersection(two)
    count1 = 0
    count2 = 0
    for i in ipr_dict1.keys():
        if i in ipr_dict2.keys():
            #print "2", ipr_dict2[i]
            count2 += 1 
        elif i not in ipr_dict2.keys():
            #print "1", ipr_dict1
            count1 +=1
        else:
            print "hello"
            
    print "\n1",count1, "\n2",count2,"\n"

    print af.name, "IPRsdict1:", len(ipr_dict1)
    print bf.name, "IPRsdict2:", len(ipr_dict2)

    
    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", af.name, bf.name
    print "1 Union 2 = ", "\t", len(one.union(two))
    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", af.name, bf.name
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - one.intersection(two))
    #**set(['d'])**

    print "\n#####################################"
    
    print "\nGOs in common: ", af.name, bf.name
    print "1 Intersection 2 = ", "\t", len(three.intersection(four)) # 10
    #print three.intersection(four)
    print "uniq ac go1:", uniq_gos1_ac[:10], "\nTOTAL UNIQUE ", len(uniq_gos1_ac)
    print "uniq ac go2:", uniq_gos2_ac[:10], "\nTOTAL UNIQUE ", len(uniq_gos1_ac)

    count_1 = 0
    count_2 = 0
    for i in go_dict1.keys():
        if i in go_dict2.keys():
            #print "2", go_dict2[i]
            count_2 += 1 
        elif i not in go_dict2.keys():

            #print "1", go_dict1
            count_1 +=1
        else:
            print "hello"
            
    print "\nonly 1:",count_1, "\nboth:",count_2
    print af.name, "GOsdict1:", len(go_dict1)
    print bf.name, "GOsdict2:", len(go_dict2)

    #**set(['a', 'c', 'b'])**
    print "\nAll GOs ", af.name, bf.name
    print "1 Union 2 = ", "\t", len(three.union(four)) # 157
    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nGOs ", af.name, bf.name
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(three.union(four)  - three.intersection(four)) #  147
    #**set(['d'])**



#get_enrichment_with_bg("2pep_260.txt", "2pep_246_down.txt") # output/5671_bmq_ac_ipr_go.txt", "output/5661_bmq_ac_ipr_go.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")

#get_enrichment_with_bg("output/venns/5664_bmq_ac_ipr_go.txt", "output/venns/5661_bmq_ac_ipr_go.txt")


#read_tab_enrichment_ipr("output/5671_bmq_ac_ipr_go.txt", "output/5661_bmq_ac_ipr_go.txt") #"outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt")


#############################################################

def set_figure(fig):
    
    import numpy as np
    from matplotlib_venn import venn3, venn3_circles

    ax=fig.add_subplot(1,1,1)
    v = venn3(subsets=(1, 1, 1, 1, 1, 1, 1), set_labels = ('A', 'B', 'C'),ax=ax)
    v.get_patch_by_id('100').set_alpha(1.0)
    v.get_patch_by_id('100').set_color('white')
    v.get_label_by_id('100').set_text('Unknown')
    v.get_label_by_id('A').set_text('Set "A"')
    c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed',ax=ax)
    c[0].set_lw(1.0)
    c[0].set_ls('dotted')
    ax.set_title("Sample Venn diagram")
    ax.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
                ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    #show()

#set_figure(fig)

###############################################################

def venn_chart3():

    # http://w3facility.org/question/python-matplotlib-venn-diagram/

    import pylab as plt
    #from matplotlib_venn import *

    v = venn3(subsets=(1,1,0,1,0,0,0))
    v.get_label_by_id('100').set_text('First')
    v.get_label_by_id('010').set_text('Second')
    v.get_label_by_id('001').set_text('Third')

    plt.title("Not a Venn diagram")
    plt.show()

#venn_chart3()

##############################################################3

def venn_csv():

    import csv

    f = open("sample-sales.csv",'rt')
    reader = csv.reader(f)
    shoes = set()
    belts = set()
    shirts = set()
    for row in reader:
        customer = (row[0],row[1])
        category = row[3]
        if category == "Shoes":
            shoes.add(customer)
        if category == "Belt":
            belts.add(customer)
        if category == "Shirt":
            shirts.add(customer)

    f.close()

    print "%s customers have purchased shoes" % len(shoes)
    print "%s customers have purchased belts" % len(belts)
    print "%s customers have purchased shoes but not belts" % len(shoes - belts)
    print "%s customers have purchased shoes and belts" % len(shoes & belts)
    print "%s customers have purchases shoes and shirts" % len(shoes & shirts)
    print "%s customers have purchased shoes, belts and shirts" % len(shoes & belts & shirts)
    print "The following customers are our most valued. They have purchased shoes & belts & shirts:"

    for customer in shoes & belts & shirts:
        print customer

#venn_csv()
###########################################################3

def venn_chart_1():#aa,ab,ac,bb,bc,cc,abc):

    # http://matthiaseisen.com/pp/patterns/p0145/g    

    # Subset sizes
    s = (
        2,    # Abc
        3,    # aBc

        4,    # ABc
        3,    # abC
        1,    # AbC
        0.5,  # aBC
        4,    # ABC
        )

    v = venn3(subsets=s, set_labels=('A', 'B', 'C'))

    # Subset labels
    v.get_label_by_id('100').set_text('Abc')
    v.get_label_by_id('010').set_text('aBc')
    v.get_label_by_id('110').set_text('ABc')
    v.get_label_by_id('001').set_text('Abc')
    v.get_label_by_id('101').set_text('aBc')
    v.get_label_by_id('011').set_text('ABc')
    v.get_label_by_id('111').set_text('ABC')

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

    plt.show()

# venn_chart_1()

##############################################################3
from pylab import *

def pie_chart1():

    # http://matplotlib.org/1.2.1/examples/pylab_examples/pie_demo.html

    """
    Make a pie chart - see
    http://matplotlib.sf.net/matplotlib.pylab.html#-pie for the docstring.

    This example shows a basic pie chart with labels optional features,
    like autolabeling the percentage, offsetting a slice with "explode",
    adding a shadow, and changing the starting angle.


    """
    

    # make a square figure and axes
    figure(1, figsize=(6,6))
    ax = axes([0.1, 0.1, 0.8, 0.8])

    # The slices will be ordered and plotted counter-clockwise.
    labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    fracs = [15, 30, 45, 10]
    explode=(0, 0.05, 0, 0)

    pie(fracs, explode=explode, labels=labels,
                        autopct='%1.1f%%', shadow=True, startangle=90)
                    # The default startangle is 0, which would start
                    # the Frogs slice on the x-axis.  With startangle=90,
                    # everything is rotated counter-clockwise by 90 degrees,
                    # so the plotting starts on the positive y-axis.


    title('Raining Hogs and Dogs', bbox={'facecolor':'0.8', 'pad':5})
    show()

#pie_chart1()
    
################################################################
######     CHARTS
###############################################################

def pie_chart_colors():

    # http://matplotlib.org/1.4.0/examples/pie_and_polar_charts/pie_demo_features.html

    """
    Demo of a basic pie chart plus a few additional features.
    In addition to the basic pie chart, this demo shows a few optional features:

        * slice labels
            * auto-labeling the percentage
                * offsetting a slice with "explode"
                    * drop-shadow
                        * custom start angle
    Note about the custom start angle:

    The default ``startangle`` is 0, which would start the "Frogs" slice on the
    positive x-axis. This example sets ``startangle = 90`` such that everything is
    rotated counter-clockwise by 90 degrees, and the frog slice starts on the
    positive y-axis.
    """

    import matplotlib.pyplot as plt

    # The slices will be ordered and plotted counter-clockwise.
    labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    sizes = [15, 30, 45, 10]
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    explode = (0, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
                    autopct='%1.1f%%', shadow=True, startangle=90)
    # Set aspect ratio to be equal so that pie is drawn as a circle.
    plt.axis('equal')

    plt.show()


# works
#"pie_chart_colors()


def pie_chart():

    names = []
    nums = []
    af = open("A0AAK3_sum_counts.txt","r")
    for i in af.readlines():
        #print i
        i = i.strip("\n\n")
        #print i
        #for i in i.split("\t"):
        line = i.split("\t")
        if len(line) == 2:
            for i in line:
                if len(i.split(" ")) > 1:
                    names.append(line[0])
                    nums.append(line[1])
                    #for i in line: print i

    print names
    print nums
    

    # The slices will be ordered and plotted counter-clockwise.
    #labels = 'Frogs', 'Hogs', 'Dogs', 'Logs'
    #l = ['a', 2, 'c']
    labels = str(names)[1:-1]
    print labels #= names
    
    #sizes = [15, 30, 45, 10]
    sizes = [int(x) for x in nums]
    print sizes #= nums
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral']
    explode = (0, 0.1, 0, 0) # only "explode" the 2nd slice (i.e. 'Hogs')

    plt.pie(sizes, explode=explode, labels=labels, colors=colors,
            autopct='%1.1f%%', shadow=True, startangle=90)

    # Set aspect ratio to be equal so that pie is drawn as a circle.

    plt.axis('equal')
    plt.show()

#################################################################

###################################################################3
#### SUMMARY MANUAL
#################################################################3

sourceFile = "ParentChildTreeFile_2908.txt"

def child2family(sourceFile, filterFile):

    from collections import defaultdict

    # From FILE finds IPR ID and adresses level #

    newFile = str(filterFile[:-4])+"domain_levels.txt"
    
    ac2ipr_dict = defaultdict(list)

    source_file = open(sourceFile, "rt")
    filter_file = open(filterFile, "rt")
    filter_list = []
    acs1 = []
    for aline in filter_file:
        aline = aline.strip().split("\t")
        ac1 = aline[0]
        for i in aline:
            if i[:3] == "IPR":
                ipr = i
                #print ipr
                ac2ipr_dict[ac1].append(ipr)
                if ac1 not in acs1:
                    acs1.append(ac1)
                #filter_list.append(i.strip())
                filter_list.append(ipr)
    new_file = open( newFile , 'w', 0)
    print >> new_file, 'AC\tDomain\tDomain Family\tLevel' # new_file.write("Domain\tDomain Family\tLevel")
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


#filterFile = "source2.txt"
#filterFile = "tax5671_bmq_ac_ipr_go_wsa.txt"


#child2family(sourceFile, filterFile, newFile)

    

time_end = datetime.datetime.fromtimestamp(time.time())
print("Time elapsed: ", str(time_end - time_begin))

