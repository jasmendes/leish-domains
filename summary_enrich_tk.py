
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

#import matplotlib.pyplot as plt

#import pylab as plt

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

from summary_plot_tk import *

import winsound

global theInputFiles

from pylab import figure, show, legend, ylabel


###################################################################################
#### PLOT   GRAPHS               SUMMARY                                       ###########
###################################################################################





def plot_2ylines(df):

    # http://thomas-cokelaer.info/blog/2012/04/481/

    from collections import OrderedDict
    from collections import defaultdict

    global theInputFiles
    #time_begin = datetime.datetime.fromtimestamp(time.time())

    

    # Open data to plot

    dfile = open(df, "r")
    first = dfile.readline()
    first = first.strip()
    fields = first.split("\t")
    print fields

    pvals6 = []
    counts6 = []

    
    
    gos = []
    counts = []
    pvals = []
    lines = []
    lines_id = {}
    count2pv_dict = defaultdict(list)
    dict_ipr2pv = defaultdict(list)
    dict_ipr2cts = defaultdict(list)
    dict_ipr2or = defaultdict(list)
    #dict_ipr2pv = {}
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
        if len(items) == 8:
            #print items
            xlabel = items[0]
            #gos.append(items[0])
            lines.append(items)
            ipr = items[0]
            pv =items[-2]
            odds = i[-1]
            #pvals.append(items[-2])
            cts = items[2]
            #print cts, pv
            pvals6.append(pv)
            counts6.append(cts)
            
           
            #print cts
            #counts.append(items[2])
            #count2pv_dict[cts].append(pv) #= pv
            count2pv_dict[pv] = cts# append(pv)   USING DICT of P_VALUES
            joint = "\t".join(items)
            
            lines_id[ipr]= joint
            #print joint
            dict_ipr2pv[ipr] = str(pv)
            dict_ipr2cts[ipr] = str(cts)
            dict_ipr2or[ipr] = str(odds)
            
        else:
            print len(items), items
        
            

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

    #pv_sorted = OrderedDict(sorted(dict_ipr2pv.itens(), key=lambda: x: x[1]))
    #it = iter(sorted(d.iteritems()))
    #for key in sorted(dict_ipr2pv):
    #    print dict_ipr2pv[key]

    from operator import itemgetter
    
    pvals2 = []
    counts2 = []
    
    lf = open(str(dfile.name[:-4])+"_order_pv.txt","w")
    header1 = "\t".join(fields)
    lf.write(header1+"\n")

    mf = open(str(dfile.name[:-4])+"_order_or.txt","w")
    mf.write(header1+"\n")

    pvals5 = []
    counts5 = []

    pvals4 = []
    counts4 = []
    or_sorted = sorted(dict_ipr2or.iteritems(), key=lambda x:x[1]) 
    for i, j in or_sorted:
        pvals5.append(j)
        counts5.append(dict_ipr2cts[i])
        #print i, dict_ipr2cts[i]


    order_or = sorted(lines, key=itemgetter(-1)) #lines.sort(key=lambda x: x[-2]) :error
    
    for i in order_or: # list inside list
        #print i
        if len(i) == 8:
            val = i[2]  #else "999"
            #print val
            new_line2 = "\t".join(i)
            mf.write(new_line2+"\n")
            if int(val) < 1000:
                pvals4.append(i[-2])
                #print i[-2], i[2]
                counts4.append(i[2])
            else:
                print "he"
                print i[2],i[-2], val
        else:
            print i
    

    
    order_by_pv = sorted(lines, key=itemgetter(-2)) #lines.sort(key=lambda x: x[-2]) :error
    print len(order_by_pv), len(lines)
    
    pvals3 = []
    counts3 = []
    #pvs_sorted = OrderedDict(sorted(dict_ipr2pv.items(), key=lambda x: x[1]))

    pvs_sorted = sorted(dict_ipr2pv.iteritems(), key=lambda x:x[1]) # CUTS
    print len(dict_ipr2pv)
    for i, j in pvs_sorted: # ipr : pvalue : counts
        #print i, j , lines_id[i] # dict_ipr2cts[i]
        pvals3.append(j)
        counts3.append(dict_ipr2cts[i])
        #print i, dict_ipr2cts[i]

        
        
    for i in order_by_pv: # list inside list
        #print i
        if len(i) == 8:
            val = i[2]  #else "999"
            #print val
            new_line2 = "\t".join(i)
            lf.write(new_line2+"\n")
            if int(val) < 1000:
                pvals2.append(i[-2])
                
                #print i[-2], i[2]
                counts2.append(i[2])
            else:
                print "he"
                print i[2],i[-2], val
        else:
            print i

    lf.close()
    mf.close()
    #print dfile.name[:-4]
    dfile.close()
    lines = []
    #global theInputFiles
    
    theInputFiles.append(lf.name) # PV - IS NOT IMPORTING FILE TO LIST
    #theInputFiles.append(mf.name) # OR
    #print theInputFiles
    #updatecombolist(theFilenames1)

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
    legend(loc="upper left")
    #legend()

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
    legend()
    #host.xaxis_date()
    #fig1.show()
    #draw()
    show()
    
    #Freq = 1000 # Set Frequency To 2500 Hertz
    #Dur = 100 # Set Duration To 1000 ms == 1 second
    #winsound.Beep(Freq,Dur)
    #print "(complete!)"
        
    #time_end = datetime.datetime.fromtimestamp(time.time())
    #print("Time elapsed: ", str(time_end - time_begin))

    
    
#df = "output/A4I8L8_entry_ipr_name_gene_go_path_local_exi_ipr_enrichments.txt"
#df1 = "output/A4I8Y2_entry_ipr_name_gene_go_path_local_exi_ipr_enrichments.txt"
#df2 = "output/A4I8Y2_entry_ipr_name_gene_go_path_local_exi_go_enrichments.txt"

#df3 = "output/E9AGK1_259_bs_uniprot_ipr_enrichments.txt"#"output/A4I1D9_248_bs_biomart_go_enrichments.txt"
#df5 = "output/E9AGK1_259_bs_uniprot_go_enrichments.txt"#"output/A4I1D9_248_bs_biomart_ipr_enrichments.txt"7
#df5 = "output/A4HRR5_259_bs_uniprot_ipr_enrichments.txt"#"output/A4I1D9_248_bs_biomart_ipr_enrichments.txt"7

#plot_2ylines(df10)
#plot_2ylines(df5)
    
# A4I8Y2_entry_ipr_name_gene_go_path_local_exi_go_enrichments
#df6 = "output/A4I8K7_259_bs_uniprot_ipr_enrichments.txt"
#df6 = "output/A4HUZ8_259_bs_uniprot"

#df6 = "output/E9AHP1_259_bs_uniprot_ipr_enrichments.txt"
#df7 = "output/E9AHP1_259_bs_uniprot_go_enrichments.txt"

# 260
# WORKS UNIPROT
#df8 = "output/2pep_260_260_match_tax5671_bs_uniprot_ipr_enrichments.txt"
#df9 = "output/2pep_260_260_match_tax5671_bs_uniprot_go_enrichments.txt"

#df14 = "output/2pep_246_down_246_match_tax5671_bs_uniprot_ipr_enrichments.txt"
#df15 = "output/2pep_246_down_246_match_tax5671_bs_uniprot_go_enrichments.txt"

# WORKS BMQ

# 246
#df10 = "output/2pep_246_down_246_match_tax5671_bmq_ac_ipr_go_ipr_enrichments.txt"
#df11 = "output/2pep_246_down_246_match_tax5671_bmq_ac_ipr_go_go_enrichments.txt"

#df12 = "output/2pep_260_260_match_tax5671_bmq_ac_ipr_go_go_enrichments.txt"
#df13 = "output/2pep_260_260_match_tax5671_bmq_ac_ipr_go_ipr_enrichments.txt"

#plot_2ylines(df8)
#plot_2ylines(df9)

#plot_2ylines(df10)
#plot_2ylines(df11)
#plot_2ylines(df12)
#plot_2ylines(df13)

#plot_2ylines(df8)
#plot_2ylines(df9)

#plot_2ylines(df15)
#plot_2ylines(df14)


#df11 = "output/2pep_246_down_246_match_tax5671_bmq_ac_ipr_go_enrichments.txt"




################################################################################


def plot_2ylines_pval(df):

    # http://thomas-cokelaer.info/blog/2012/04/481/

    from collections import OrderedDict
    from collections import defaultdict
    import numpy as np
    #import matplotlib as plt
    import matplotlib.pyplot as plt
    

    from operator import itemgetter

    global theInputFiles
    #time_begin = datetime.datetime.fromtimestamp(time.time())

    

    # Open data to plot

    dfile = open(df, "r")
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
    #dict_ipr2pv = {}
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
        if len(items) == 8:
            #print items
            xlabel = items[0]
            #gos.append(items[0])
            lines.append(items)
            ipr = items[0]
            pv =items[ind_pv]
            odds = i[ind_or]
            
            #pvals.append(items[-2])
            cts = items[2]
            #print cts, pv
            pvals6.append(pv)
            counts6.append(cts)
            
           
            #print cts
            #counts.append(items[2])
            #count2pv_dict[cts].append(pv) #= pv
            count2pv_dict[pv] = cts# append(pv)   USING DICT of P_VALUES
            joint = "\t".join(items)
            
            lines_id[ipr]= joint
            #print joint
            dict_ipr2pv[ipr] = str(pv)
            dict_ipr2cts[ipr] = str(cts)
            dict_ipr2or[ipr] = str(odds)
            
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
    
    lf = open(str(dfile.name[:-4])+"_order_pv.txt","w")
    header1 = "\t".join(fields)
    lf.write(header1+"\n")

    mf = open(str(dfile.name[:-4])+"_order_or.txt","w")
    mf.write(header1+"\n")

    pvals5 = []
    counts5 = []

    pvals4 = []
    counts4 = []
    or_sorted = sorted(dict_ipr2or.iteritems(), key=lambda x:x[1]) 
    for i, j in or_sorted:
        pvals5.append(j)
        counts5.append(dict_ipr2cts[i])
        #print i, dict_ipr2cts[i]


    order_or = sorted(lines, key=itemgetter(ind_or)) #lines.sort(key=lambda x: x[-2]) :error
    
    for i in order_or: # list inside list
        #print i
        if len(i) == 8:
            val = i[2]  #else "999"
            #print val
            new_line2 = "\t".join(i)
            mf.write(new_line2+"\n")
            if int(val) < 2000:
                pvals4.append(i[ind_pv])
                #print i[-2], i[2]
                counts4.append(i[2])
            else:
                print "he"
                print i[2],i[ind_pv], val
        else:
            print i
    

    
    order_by_pv = sorted(lines, key=itemgetter(ind_pv)) #lines.sort(key=lambda x: x[-2]) :error
    print len(order_by_pv), len(lines)
    
    pvals3 = []
    counts3 = []
    #pvs_sorted = OrderedDict(sorted(dict_ipr2pv.items(), key=lambda x: x[1]))

    pvs_sorted = sorted(dict_ipr2pv.iteritems(), key=lambda x:x[1]) # CUTS
    print len(dict_ipr2pv)
    for i, j in pvs_sorted: # ipr : pvalue : counts
        #print i, j , lines_id[i] # dict_ipr2cts[i]
        pvals3.append(j)
        counts3.append(dict_ipr2cts[i])
        #print i, dict_ipr2cts[i]

        
    # CORRECT VERSION
    bjmn2 = []
    odds2 = []
    bfrn2 = []
    for i in order_by_pv: # list inside list
        #print i
        if len(i) == 8:
            val = i[2]  #else "999"
            #print val
            new_line2 = "\t".join(i)
            lf.write(new_line2+"\n")
            if int(val) < 1000:
                pvals2.append(i[ind_pv])
                counts2.append(i[2])
                odds2.append(i[ind_or])
                
                
                #print i[-2], i[2]
                
                
            else:
                print "hello i[2],i[ind_pv], val "
                print i[2],i[ind_pv], val
        else:
            print i

    lf.close()
    mf.close()
    #print dfile.name[:-4]
    dfile.close()
    lines = []
    #global theInputFiles
    
    theInputFiles.append(lf.name) # PV - IS NOT IMPORTING FILE TO LIST
    #theInputFiles.append(mf.name) # OR
    #print theInputFiles
    #updatecombolist(theFilenames1)



    title = fields[0] +" vs. " +fields[ind_pv] + "\n" + dfile.name
    # Create the general figure
    #fig1 = figure()
    
    fig1 = plt.figure()
    fig1.canvas.set_window_title(str(lf.name[:-4]))

   
    fig1.suptitle(str(title))

    #t = np.arange(0.0, 1.0, 0.01)
    #s = np.sin(2*np.pi*t)
    #line, = ax.plot(t, s, color='blue', lw=2)
    
    

    # and the first axes using subplot populated with data
    ax1 = fig1.add_subplot(111)
    line1 = ax1.plot(pvals2, 'o-', label='p-value')

    #line11 = ax1.plot(odds2, 'o-', label='odds ratio')
    
    #ax1.yaxis.tick_right()
    #ax1.yaxis.set_label_position("right")
    ylabel("P-value")
    legend(loc="upper left")
    #legend()

    # draw vertical / horizontal line from (70,100) to (70, 250)
    ax1.plot([0, len(counts2)], [0.05, 0.05], 'g-', lw=2)

    #line0 = ax1.plot(t, bline, color='blue', lw=2)

    # now, the second axes that shares the x-axis with the ax1

    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
    line2 = ax2.plot(counts2, 'xr-', label='counts')
    
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ylabel(" Counts")
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
    legend()
    #legend((line1, line2), ("p-value", "counts"))
    #host.xaxis_date()
    #fig1.show()
    #draw()
    show()
    
    #Freq = 1000 # Set Frequency To 2500 Hertz
    #Dur = 100 # Set Duration To 1000 ms == 1 second
    #winsound.Beep(Freq,Dur)
    #print "(complete!)"
        
    #time_end = datetime.datetime.fromtimestamp(time.time())
    #print("Time elapsed: ", str(time_end - time_begin))

    
    
# 260
# WORKS UNIPROT


#df14 = "output/2pep_246_down_246_match_tax5671_bs_uniprot_ipr_enrichments.txt"
#df15 = "output/2pep_246_down_filter_246_tax5671_bs_uniprot_ipr_enrichments.txt"

# df11 = "output/2pep_246_down_246_match_tax5671_bmq_ac_ipr_go_enrichments.txt"


#plot_2ylines(df14)
#df16 = "20150827015604\\2peps_log2dot_filter_2293_A4HRT5_filter_2294_tax5671_bs_uniprot_m_ipr_enrichments_order_pv.txt"

#plot_2ylines_pval(df16)

#####################################################################


from pylab import figure, show, legend, ylabel
import winsound



def plot_2ylines_or(df):

    # http://thomas-cokelaer.info/blog/2012/04/481/

    from collections import OrderedDict
    from collections import defaultdict
    import numpy as np
    #import matplotlib as plt
    import matplotlib.pyplot as plt
    

    from operator import itemgetter

    global theInputFiles
    #time_begin = datetime.datetime.fromtimestamp(time.time())

    

    # Open data to plot

    dfile = open(df, "r")
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
    #dict_ipr2pv = {}
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
        if len(items) == 8:
            #print items
            xlabel = items[0]
            #gos.append(items[0])
            lines.append(items)
            ipr = items[0]
            pv =items[ind_pv]
            odds = i[ind_or]
            
            #pvals.append(items[-2])
            cts = items[2]
            #print cts, pv
            pvals6.append(pv)
            counts6.append(cts)
            
           
            #print cts
            #counts.append(items[2])
            #count2pv_dict[cts].append(pv) #= pv
            count2pv_dict[pv] = cts# append(pv)   USING DICT of P_VALUES
            joint = "\t".join(items)
            
            lines_id[ipr]= joint
            #print joint
            dict_ipr2pv[ipr] = str(pv)
            dict_ipr2cts[ipr] = str(cts)
            dict_ipr2or[ipr] = str(odds)
            
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
    
    lf = open(str(dfile.name[:-4])+"_order_pv.txt","w")
    header1 = "\t".join(fields)
    lf.write(header1+"\n")

    mf = open(str(dfile.name[:-4])+"_order_or.txt","w")
    mf.write(header1+"\n")

    pvals5 = []
    counts5 = []

    pvals4 = []
    counts4 = []
    or_sorted = sorted(dict_ipr2or.iteritems(), key=lambda x:x[1]) 
    for i, j in or_sorted:
        pvals5.append(j)
        counts5.append(dict_ipr2cts[i])
        #print i, dict_ipr2cts[i]


    order_or = sorted(lines, key=itemgetter(ind_or)) #lines.sort(key=lambda x: x[-2]) :error
    
    for i in order_or: # list inside list
        #print i
        if len(i) == 8:
            val = i[2]  #else "999"
            #print val
            new_line2 = "\t".join(i)
            mf.write(new_line2+"\n")
            if int(val) < 2000:
                pvals4.append(i[ind_pv])
                #print i[-2], i[2]
                counts4.append(i[2])
            else:
                print "he"
                print i[2],i[ind_pv], val
        else:
            print i
    

    
    order_by_pv = sorted(lines, key=itemgetter(ind_pv)) #lines.sort(key=lambda x: x[-2]) :error
    print len(order_by_pv), len(lines)
    
    pvals3 = []
    counts3 = []
    #pvs_sorted = OrderedDict(sorted(dict_ipr2pv.items(), key=lambda x: x[1]))

    pvs_sorted = sorted(dict_ipr2pv.iteritems(), key=lambda x:x[1]) # CUTS
    print len(dict_ipr2pv)
    for i, j in pvs_sorted: # ipr : pvalue : counts
        #print i, j , lines_id[i] # dict_ipr2cts[i]
        pvals3.append(j)
        counts3.append(dict_ipr2cts[i])
        #print i, dict_ipr2cts[i]

        
    # CORRECT VERSION
    bjmn2 = []
    odds2 = []
    bfrn2 = []
    for i in order_by_pv: # list inside list
        #print i
        if len(i) == 8:
            val = i[2]  #else "999"
            #print val
            new_line2 = "\t".join(i)
            lf.write(new_line2+"\n")
            if int(val) < 1000:
                pvals2.append(i[ind_pv])
                counts2.append(i[2])
                odds2.append(i[ind_or])
                
                
                #print i[-2], i[2]
                
                
            else:
                print "hello i[2],i[ind_pv], val "
                print i[2],i[ind_pv], val
        else:
            print i

    lf.close()
    mf.close()
    #print dfile.name[:-4]
    dfile.close()
    lines = []
    #global theInputFiles
    
    theInputFiles.append(lf.name) #pv
    theInputFiles.append(mf.name) # OR
    #print theInputFiles
    #updatecombolist(theFilenames1)



    title = fields[0] +" vs. " +fields[ind_pv] + "\n" + dfile.name
    # Create the general figure
    #fig1 = figure()
    
    fig1 = plt.figure()
    fig1.canvas.set_window_title(str(mf.name[:-4]))

   
    fig1.suptitle(str(title))

    #t = np.arange(0.0, 1.0, 0.01)
    #s = np.sin(2*np.pi*t)
    #line, = ax.plot(t, s, color='blue', lw=2)
    
    

    # and the first axes using subplot populated with data
    ax1 = fig1.add_subplot(111)
    line1 = ax1.plot(pvals2, 'o-', label='p-value')

    line11 = ax1.plot(odds2, 'o-', label='odds ratio')
    
    #ax1.yaxis.tick_right()
    #ax1.yaxis.set_label_position("right")
    ylabel("P-value | Odds Ratio")
    legend(loc="upper left")
    #legend()

    # draw vertical line from (70,100) to (70, 250)
    ax1.plot([0, len(counts2)], [0.05, 0.05], 'g-', lw=2)

    #line0 = ax1.plot(t, bline, color='blue', lw=2)

    # now, the second axes that shares the x-axis with the ax1

    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
    line2 = ax2.plot(counts2, 'xr-', label='counts')
    
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ylabel("Counts")
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
    legend()
    #legend((line1, line2), ("p-value", "counts"))
    #host.xaxis_date()
    #fig1.show()
    #draw()
    show()
    
    #Freq = 1000 # Set Frequency To 2500 Hertz
    #Dur = 100 # Set Duration To 1000 ms == 1 second
    #winsound.Beep(Freq,Dur)
    #print "(complete!)"
        
    #time_end = datetime.datetime.fromtimestamp(time.time())
    #print("Time elapsed: ", str(time_end - time_begin))

    
    
# 260
# WORKS UNIPROT


#df14 = "output/2pep_246_down_246_match_tax5671_bs_uniprot_ipr_enrichments.txt"
#df15 = "output/2pep_246_down_filter_246_tax5671_bs_uniprot_ipr_enrichments.txt"

# df11 = "output/2pep_246_down_246_match_tax5671_bmq_ac_ipr_go_enrichments.txt"


#plot_2ylines(df14)
#df16 = "20150827015604\\2peps_log2dot_filter_2293_A4HRT5_filter_2294_tax5671_bs_uniprot_m_ipr_enrichments_order_pv.txt"

#plot_2ylines_or(df16)

#################################################################################

####################################################################################3
#####  ENRICHMENT ANALYSIS - IPR AND GO -
#####################################################################################3

import scipy as sp
import scipy.stats
import numpy as np


def enrichment_analysis_paths(af, bf):

    print "\nWork in process...."

    """Path Enrichment Analysis                                    """


   
    global theInputFiles
    
    print "choice 1 = ", af
    print "choice 2 = ", bf
    
    cfile = str(af[:-4])+"_enrich_path.txt"
    cf = open(cfile, "w")
    
    
    iprs1 = []
    iprs2 = []
    gos1 = []
    gos2 = []
    iprn_dict = {}
    gon_dict = {}
    path_dict = {}
    tot_paths = []
    all_paths = []




    
    
    acs1 = []

    iacs1 = []
    iprc1 = 0
    goc1 = 0
    iprs11 = []
    gos11 = []
    all_ac1 = []
    all_ac2 = []

    pacs1 = []
    path_acs1 = []
    pacs11 = []
    import collections
    
    num_lines1 = sum(1 for line in open(af))
    with open(af,"r") as f1:
        for i in f1.readlines():
            #print i
            i = i.strip()
            items = i.split("\t")
            ac = items[0]
            if ac not in all_ac1:
                all_ac1.append(ac)
            for j in items:
                #print j
                if j[0:8] in all_paths:
                    path = j.split(":")[0]
                    pathname = j.split(":")[1]
                    if ac not in path_acs1:
                        path_acs1.append(ac)
                    pacs1.append(path)
                    if path not in pacs11:
                        pacs11.append(path)
                
                

    pacs2 = []
    pacs22 = []
    path_acs2 = []
    path_ac_dict =  defaultdict(list)
    path_dictn = []
    
    with open(bf,"r") as f2:
        for i in f2.readlines():
            #print i
            i = i.strip()
            items = i.split("\t")
            ac = items[0]
            if ac not in all_ac2:
                all_ac2.append(ac)
            for j in items:
                #print j
                if j[0:8] in all_paths:
                    path_dictn[path] = pathname
                    path_ac_dict[ac] = j
                    path = j.split(":")[0]
                    pathname = j.split(":")[1]
                    if ac not in path_acs2:
                        path_acs2.append(ac)
                    pacs2.append(path)
                    if path not in pacs22:
                        pacs22.append(path)
                
                

    

    ##############3   VENN DIAGGRAM   IPR1 vs. IPR2 #######################3

    print "TOTAL Paths1 : ",len(pacs1),"\nTOTAL Paths2 : " , len(pacs2)
    

    print "UNIQUE Paths1 : ",len(pacs11),"\nUNIQUE Paths 2 : " , len(pacs22)
    


    ################################    VENN    ###################3
    one = set(pacs1) #uniq_iprs1)
    two = set(pacs2) # uniq_iprs2)

    # IPR IDs
    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", af, bf
    print "1 Union 2 = ", "\t", len(one.union(two))
    #a = most.common(one)
    
    print "\nDomains in common: ", af, bf
    print "1 Intersection 2 = ", "\t", len(one.intersection(two))
    #print one.intersection(two)

    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", af, bf
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - one.intersection(two))
    #**set(['d'])**

    #a = [1, 2, 3, 4, 5]
    #b = [9, 8, 7, 6, 5]
    #i for i, j in zip(a, b) if i == j]

    
    
    ############ IPR SUM ##########
    
    print"\nTOTAL ACS 1 :",len(all_ac1),
    print"\tTOTAL ACS 2 :",len(all_ac2)
    
    c  = Counter(pacs1)
    print "\nTOTAL iprs1 found : ", len(pacs1) # 6597
    print "UNIQUE iprs1 : ", len(c), "\n", c.most_common(10) # 2088
    #for i in dict(c.most_common()):
    #   cf.write(str(i)+"\t"+str(c[i])+"\n")
    # print i, c[i]

    c_dict = dict(c.most_common())
        
    ### COVERAGE IPR 1
    #print iprc1, acs1
    try:
        cov =  (len(pacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enr = (len(pacs11)/float(len(pacs1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation IPR 1 :\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No IPRs 1 found"
    

    
    d  = Counter(pacs2)
    print "\nTOTAL IPRs2 found : ", len(pacs2) # 6597
    print "UNIQUE IPRs2 : ", len(d), "\n", d.most_common(10) # 3606
    #for i in dict(d.most_common()):
       # cf.write(str(i)+"\t"+str(d[i])+"\n")
        #print i, d[i]
    d_dict = dict(d.most_common())

    ### COVERAGE IPR 2
    #print iprc2, acs2
    try:
        cov = (len(pacs2)/float(len(all_ac2))) # most common(uncharacterized) / total found
        enr = (len(pacs22)/float(len(pacs2))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation IPR 2 :\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No IPRs 2 found!"
        #e  = Counter(ippr_names)
    #print "\nTOTAL IPRs NAMES found : ", len(ippr_names) # 6597
    #print "UNIQUE IPR NAMES : ", len(e), "\n", e.most_common(10) # 3606
    #for i in dict(e.most_common(10)):
     #   print i, e[i]


    

    #for i in g_dict:
        #cf.write(str(i)+"\t"+str(g[i])+"\n")
        #print i, g_dict[i]
        #oddsratio, pvalue = stats.fisher_exact([[1100,6848], [11860,75292]])

    ### PVALUE + ODDSRATIO
    
    #print iprn_dict["IPR030470"]
    #iprn_dict["IPR030470"] = "null"
    df.write("Path ID\tPath Name\tcounts\ttotal_counts\tquery\tall\tp-value\todds_ratio\n") # lack id name. how to? fromgooripr_insert bmq_ipr2go
    results = []
    ipr_pv_dict = {}
    pvalues_ipr = []
    rank = 0
    for j in sorted(dict(c.most_common()), key=lambda x :dict(c.most_common())[x], reverse=True):
        #print j, f[j]
        if j in d_dict.keys():
            #print j, g_dict[j]
            n1 = c[j]
            n2 = d_dict[j]
            #n3 = len(c)- n1
            n3 = len(pacs1) - n1  # n4 = len(d_dict) - n2
            
            n4 = len(pacs2) - n2  # print j, n1, n2, n3, n4
            #n3 = len(iacs1)
            #n4 = len(iacs2)
            
            oddsratio, pvalue = stats.fisher_exact([[n1,n2], [n3,n4]],alternative="greater")
            #print oddsratio, pvalue
            #print '%f  ' % oddsratio
            #print '%f ' % pvalue

            pvalue = "{0:.10f}".format(pvalue)
            oddsratio = "{0:.10f}".format(oddsratio)
            #print pvalue, oddsratio
            #oddsratio = remove_exponent(oddsratio)
            #pvalue = remove_exponent(pvalue)
            #print oddsratio, pvalue
            #if n1 and n2 and n3 and n4 > 0 else "null"

            #import numpy as np
            #import scipy.stats as stats

            #G = stats.poisson(30)
            #print(G.dist.expect(lambda x: (x+1), G.args, lb=0, ub=np.inf))
            # 31.0
            
            if j == "":
                pass
            if j not in iprn_dict.keys():
                iprn_dict[j] = "null"
            pvalues_ipr.append(pvalue)
            ipr_pv_dict[j]=pvalue
            
            #pvalue = lambda pvalue: float(pvalue.replace('.',','))
            #oddsratio = oddsratio = lambda oddsratio: float(oddsratio.replace('.',',')) 
            
            #print j, n1 , n2, n3, n4,"\tp-value = ", pvalue, "\tOdds Ratio = ", oddsratio
            #pvalue = float(pvalue)
            #oddsratio = float(oddsratio)
            
            
            #pvalue = abs(pvalue)
            #oddsratio = abs(oddsratio)
            #oddsratio = 1/float(oddsratio)
            
            results.append(str(pvalue)+"\t"+str(oddsratio))
            df.write(str(j)+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)
                     +"\t"+str(pvalue)+"\t"+str(oddsratio)+"\n")
            

    
    
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


    cf.close()
    df.close()


    # Plot 2y Scatter plot : ID vs p-value
    # Some points are missing while running the plot_2ylines riht away, so
    # we will split it in another function -> Scatter Enrich

    #print cfile, dfile
    #plot_2ylines2(cfile)
    #plot_2ylines2(dfile)
    #global theInputFiles
    theInputFiles.append(cf.name)
    
    #theFilenames1.append(cf.name)
    #global theFilenames1
    #theFilenames1.append(cf.name)
    #theFilenames1.append(df.name)
    
    #updatecombolist2(theInputFiles)
    

    # False Discovery Rate of the p-value 
        

    print cf.name, " saved!"
    
    #print ef, " saved!"
    #ef.close()
    #cff.close()
    #dff.close()
    
    #name1 = os.path.basename(cf)
    #name2 = os.path.basename(df)

    #global theInputFiles
    #theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)
    #Freq = 1000 # Set Frequency To 2500 Hertz
    #Dur = 100 # Set Duration To 1000 ms == 1 second
    #winsound.Beep(Freq,Dur)
    print "(complete!)"

    
#af = "output/tax5671_bs_uniprot_filter_892_A4HZW1_filter_893_kegg_lin.txt"
#bf = "_experimental/all_FC_filter_5213_kegg_lin.txt"
#enrichment_analysis_paths(af, bf)

# 347515 L major
def enrichment_analysis_domains(af, bf, kf): # n = 1000

    """Using the most common dict (COLLECTIONS) file from last functions,
    we are able to build a default proteome set (bf) 
    and compare it with input proteome set (af) 
    combining de ipr2go EBI dictionary (kf)                                     """


    # Gets tax id and retrieves uni ids list from uniprot web
    # then enrichment analysis with list of accessions

    #from os.path import basename
    # now you can call it directly with basename
    #print basename(af)
    #str(af[:-4])+
    global theInputFiles
    
    print "choice 1 = ", af
    print "choice 2 = ", bf
    print "dict ipr2go = " , kf
    cfile = str(af[:-4])+"_enrich_gos.txt"
    cf = open(cfile, "w")
    dfile = str(af[:-4])+"_enrich_iprs.txt"
    df = open(dfile, "w")
    
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
    num_lines1 = sum(1 for line in open(af))
    with open(af,"r") as f1:
        for i in f1.readlines():
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
                
                

    iacs2 = []
    #acs2 = []
    iprc2 = 0
    goc2 = 0
    iprs22 = []
    gos22 = []
    #num_lines2 = sum(1 for line in open(bf))
    with open(bf,"r") as f2:
        for i in f2.readlines():
            #print i
            i = i.strip()
            items = i.split("\t")
            ac = items[0]
            if ac not in all_ac2:
                all_ac2.append(ac)
            for j in items:
                #print j
                if j[0:3] == "IPR":
                    iprc2 += 1
                    if ac not in iacs2:
                        iacs2.append(ac)
                    #print j.split("; ")
                    for i in j.split("; "):
                        iprs2.append(i)
                        if i not in iprs22:
                            iprs22.append(i)
                elif j[0:3] == "GO:":
                    goc2 += 1
                    #print j.split("; ")
                    if ac not in gacs2:
                        gacs2.append(ac)
                    for i in j.split("; "):
                        gos2.append(i)
                        if i not in gos22:
                            gos22.append(i)
                
                

    

    ##############3   VENN DIAGGRAM   IPR1 vs. IPR2 #######################3

    print "TOTAL IPRs1 : ",len(iprs1),"\nTOTAL IPRs2 : " , len(iprs2)
    print "TOTAL GOs1 : ",len(gos1),"\nTOTAL GOs2 : " , len(gos2)

    print "UNIQUE IPRs1 : ",len(iprs11),"\nUNIQUE IPRs2 : " , len(iprs22)
    print "UNIQUE GOs1 : ",len(gos11),"\nUNIQUE GOs2 : " , len(gos22)


    ################################    VENN    ###################3
    one = set(iprs1) #uniq_iprs1)
    two = set(iprs2) # uniq_iprs2)

    # IPR IDs
    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", af, bf
    print "1 Union 2 = ", "\t", len(one.union(two))
    #a = most.common(one)
    
    print "\nDomains in common: ", af, bf
    print "1 Intersection 2 = ", "\t", len(one.intersection(two))
    #print one.intersection(two)

    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", af, bf
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - one.intersection(two))
    #**set(['d'])**

    #a = [1, 2, 3, 4, 5]
    #b = [9, 8, 7, 6, 5]
    #i for i, j in zip(a, b) if i == j]

    
    # GO TERMS
    three = set(gos1) #uniq_iprs1)
    four = set(gos2) # uniq_iprs2)

    #**set(['a', 'c', 'b'])**
    print "\nAll GOs ", af, bf
    print "2 Union 1 = ", "\t", len(four.union(three))
    #print four.union(three)
    
    print "\nGOs in common: ", af, bf
    print "2 Intersection 1 = ", "\t", len(four.intersection(three))
    #print one.intersection(two)

    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nGOs ", af, bf
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(four.union(three)  - four.intersection(three))

    ############ IPR SUM ##########
    
    print"\nTOTAL ACS 1 :",len(all_ac1),
    print"\tTOTAL ACS 2 :",len(all_ac2)
    
    c  = Counter(iprs1)
    print "\nTOTAL iprs1 found : ", len(iprs1) # 6597
    print "UNIQUE iprs1 : ", len(c), "\n", c.most_common(10) # 2088
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
    

    
    d  = Counter(iprs2)
    print "\nTOTAL IPRs2 found : ", len(iprs2) # 6597
    print "UNIQUE IPRs2 : ", len(d), "\n", d.most_common(10) # 3606
    #for i in dict(d.most_common()):
       # cf.write(str(i)+"\t"+str(d[i])+"\n")
        #print i, d[i]
    d_dict = dict(d.most_common())

    ### COVERAGE IPR 2
    #print iprc2, acs2
    try:
        cov = (len(iacs2)/float(len(all_ac2))) # most common(uncharacterized) / total found
        enr = (len(iprs22)/float(len(iprs2))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation IPR 2 :\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No IPRs 2 found!"
        #e  = Counter(ippr_names)
    #print "\nTOTAL IPRs NAMES found : ", len(ippr_names) # 6597
    #print "UNIQUE IPR NAMES : ", len(e), "\n", e.most_common(10) # 3606
    #for i in dict(e.most_common(10)):
     #   print i, e[i]


    ### GO SUM ###
    ### COVERAGE GO 1 
    
    f  = Counter(gos1)
    print "\nTOTAL GO IDS 1 found : ", len(gos1) # 6597
    print "UNIQUE GO IDs 1: ", len(f), "\n", f.most_common(10) # 1183
    #for i in dict(f.most_common()):
        #cf.write(str(i)+"\t"+str(f[i])+"\n")
        # print i, f[i]

    try:
        cov = (len(gacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enr = (len(gos11)/float(len(gos1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 1 :\n\tAc Coverage : ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No GOs 1 found!"
    

    g  = Counter(gos2)
    print "\nTOTAL GO 2 found : ", len(gos2)# 6597
    print "UNIQUE GO 2: ", len(g), "\n", g.most_common(10) # 1183
    g_dict = dict(g.most_common())

    
    ### COVERAGE GO 2
    try:
        cov = (len(gacs2)/float(len(all_ac2))) # most common(uncharacterized) / total found
        enr = (len(gos22)/float(len(gos2))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 2:\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No GOs 2 found!"



    #for i in g_dict:
        #cf.write(str(i)+"\t"+str(g[i])+"\n")
        #print i, g_dict[i]
        #oddsratio, pvalue = stats.fisher_exact([[1100,6848], [11860,75292]])

    ### PVALUE + ODDSRATIO
    
    #print iprn_dict["IPR030470"]
    #iprn_dict["IPR030470"] = "null"
    df.write("IPR ID\tIPR Name\tcounts\ttotal_counts\tquery\tall\tp-value\todds_ratio\n") # lack id name. how to? fromgooripr_insert bmq_ipr2go
    results = []
    ipr_pv_dict = {}
    pvalues_ipr = []
    rank = 0
    for j in sorted(dict(c.most_common()), key=lambda x :dict(c.most_common())[x], reverse=True):
        #print j, f[j]
        if j in d_dict.keys():
            #print j, g_dict[j]
            n1 = c[j]
            n2 = d_dict[j]
            #n3 = len(c)- n1
            n3 = len(iprs1) - n1  # n4 = len(d_dict) - n2
            
            n4 = len(iprs2) - n2  # print j, n1, n2, n3, n4
            #n3 = len(iacs1)
            #n4 = len(iacs2)
            
            oddsratio, pvalue = stats.fisher_exact([[n1,n2], [n3,n4]],alternative="greater") #  alternative="two-sided"); "less"
            #print oddsratio, pvalue
            #print '%f  ' % oddsratio
            #print '%f ' % pvalue

            pvalue = "{0:.10f}".format(pvalue)
            oddsratio = "{0:.10f}".format(oddsratio)
            #print pvalue, oddsratio
            #oddsratio = remove_exponent(oddsratio)
            #pvalue = remove_exponent(pvalue)
            #print oddsratio, pvalue
            #if n1 and n2 and n3 and n4 > 0 else "null"

            #import numpy as np
            #import scipy.stats as stats

            #G = stats.poisson(30)
            #print(G.dist.expect(lambda x: (x+1), G.args, lb=0, ub=np.inf))
            # 31.0
            
            if j == "":
                pass
            if j not in iprn_dict.keys():
                iprn_dict[j] = "null"
            pvalues_ipr.append(pvalue)
            ipr_pv_dict[j]=pvalue
            
            #pvalue = lambda pvalue: float(pvalue.replace('.',','))
            #oddsratio = oddsratio = lambda oddsratio: float(oddsratio.replace('.',',')) 
            
            #print j, n1 , n2, n3, n4,"\tp-value = ", pvalue, "\tOdds Ratio = ", oddsratio
            #pvalue = float(pvalue)
            #oddsratio = float(oddsratio)
            
            
            #pvalue = abs(pvalue)
            #oddsratio = abs(oddsratio)
            #oddsratio = 1/float(oddsratio)
            
            results.append(str(pvalue)+"\t"+str(oddsratio))
            df.write(str(j)+"\t"+str(iprn_dict[j])+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)
                     +"\t"+str(pvalue)+"\t"+str(oddsratio)+"\n")
            

    cf.write("GO ID\tGO Term\tcounts\ttotal_counts\tquery\tall\tp-value\todds_ratio\n") # lack id name. how to? fromgooripr_insert bmq_ipr2go
    results = []
    go_pv_dict = {}
    pvalues_go = []
    for j in sorted(dict(f.most_common()), key=lambda x :dict(f.most_common())[x], reverse=True):
        #print j, f[j]
        if j in g_dict.keys():
            #print j, g_dict[j]
            n1 = f[j]
            n2 = g_dict[j]
            n3 = len(gos1) - n1  # all gos
            n4 = len(gos2) - n2 #
            #n3 = len(gacs1) - n1  # all gos
            #n4 = len(gacs2) - n2 # 
            
            #print j,f[j], g_dict[j], len(f), len(g_dict)
            #print n1, n2, n3, n4
            oddsratio, pvalue = stats.fisher_exact([[n1,n2], [n3,n4]],alternative="greater")
            #print float(oddsratio), float(pvalue)

            # "%.16f" % f if f >= 1e-16 else "0.0"
            pvalue = "{0:.10f}".format(pvalue)
            oddsratio = "{0:.10f}".format(oddsratio)
            #result_arr = fisherExact(np.array([[5, 0], [1, 4]]))
            #print result_arr
            
            # SCIPY pvalue, uses hypergeometric distribution (not efficient)
            #x = [[n1,n2],[n3,n4]]
            #result_pv = sp.stats.fisher_exact(x)
            # print result_pv[1], "pv"
            #pvalue2 = result_pv[1] # for i in result_pv:
            #odds2 = result_pv[0]
            
            #print float(oddsratio), float(pvalue)


            if j == "":
                pass
            if j not in gon_dict.keys():
                gon_dict[j] = "null"
                #break

            #pval_go = myformat(pvalue)
            pvalues_go.append(pvalue)
            go_pv_dict[j] = pvalue

            #pvalue = lambda pvalue: float(pvalue.replace('.',',')) # convertfloat = lambda x: float(str(x).replace(',','.')) 
            #oddsratio = lambda oddsratio: float(oddsratio.replace('.',',')) # convertfloat = lambda x: float(str(x).replace(',','.')) 
            
            #print j, n1 , n2, n3, n4,"\tp-value = ", pvalue, "\tOdds Ratio = ", oddsratio

            
            #pvalue = abs(pvalue)
            #oddsratio = abs(oddsratio)
            #oddsratio = 1/float(oddsratio)
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


    cf.close()
    df.close()


    # Plot 2y Scatter plot : ID vs p-value
    # Some points are missing while running the plot_2ylines riht away, so
    # we will split it in another function -> Scatter Enrich

    #print cfile, dfile
    #plot_2ylines2(cfile)
    #plot_2ylines2(dfile)
    #global theInputFiles
    theInputFiles.append(cf.name)
    theInputFiles.append(df.name) # NOT WORKING THE UPDATE
    #theFilenames1.append(cf.name)
    #global theFilenames1
    #theFilenames1.append(cf.name)
    #theFilenames1.append(df.name)
    
    #updatecombolist2(theInputFiles)
    

    # False Discovery Rate of the p-value 
        

    print cf.name, " saved!"
    print df.name, " saved!"
    #print ef, " saved!"
    #ef.close()
    #cff.close()
    #dff.close()
    kfile.close()
    #name1 = os.path.basename(cf)
    #name2 = os.path.basename(df)

    #global theInputFiles
    #theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)
    #Freq = 1000 # Set Frequency To 2500 Hertz
    #Dur = 100 # Set Duration To 1000 ms == 1 second
    #winsound.Beep(Freq,Dur)
    print "(complete!)"
    
#enrichment_analysis_domains("2pep_246_down_246_match_tax5671_bmq_ac_ipr_go.txt","output/tax5671_bmq_ac_ipr_go.txt","output/_bmq_ipr_go.txt")
#2pep_260_260_match_tax5671_bmq_ac_ipr_go.txt","output/tax5671_bmq_ac_ipr_go.txt","output/_bmq_ipr_go.txt")

#A4I5Q3_3990_bs_uniprot
# 246 vs all
#output/A4IBW8_entry_ipr_name_gene_go_path_local_exi.txt output/A4HU44_entry_ipr_name_gene_go_path_local_exi.txt

#enrichment_analysis_domains("output/2pep_246_down_filter_246_tax5671_bs_uniprot.txt","output/tax5671_bs_uniprot.txt", "output/_bmq_ipr_go.txt")



#enrichment_analysis_domains("2pep_246_down_246_match_tax5671_bmq_ac_ipr_go.txt","2peps4000_3998_match_tax5671_bmq_ac_ipr_go.txt","output/_bmq_ipr_go.txt")
#enrichment_analysis_domains("output/tax5671_bs_uniprot_filter_400_E9AG94.txt" , "output/tax5671_bs_uniprot_filter_5200_A4I2Z3.txt", "output/_bmq_ipr_go.txt")
"""
output/A4IBW8_entry_ipr_name_gene_go_path_local_exi.txt  New Selection 1 !
choice 1 =  output/A4IBW8_entry_ipr_name_gene_go_path_local_exi.txt
choice 2 =  output/A4I2I6_entry_ipr_name_gene_go_path_local_exi.txt
dict ipr2go =  output/_bmq_ipr_go.txt
TOTAL IPRs1 :  610 
TOTAL IPRs2 :  11239
TOTAL GOs1 :  252 
TOTAL GOs2 :  4971"""


###################################################################################


def enrichment_analysis_proteins(af, bf, kf): # n = 1000

    """Using the most common dict (COLLECTIONS) file from last functions,
    we are able to build a default proteome set (bf) 
    and compare it with input proteome set (af) 
    combining de ipr2go EBI dictionary (kf)                                     """


    # Gets tax id and retrieves uni ids list from uniprot web
    # then enrichment analysis with list of accessions
    print "################   Enrichment Analysis   #################"
    print "This method uses as  backgroud the total of proteins with at least one annotation. "

    #from os.path import basename
    # now you can call it directly with basename
    #print basename(af)
    #str(af[:-4])+
    global theInputFiles
    
    print "choice 1 = ", af
    print "choice 2 = ", bf
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
    num_lines1 = sum(1 for line in open(af))
    with open(af,"r") as f1:
        for i in f1.readlines():
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
                
                

    iacs2 = []
    #acs2 = []
    iprc2 = 0
    goc2 = 0
    iprs22 = []
    gos22 = []
    #num_lines2 = sum(1 for line in open(bf))
    with open(bf,"r") as f2:
        for i in f2.readlines():
            #print i
            i = i.strip()
            items = i.split("\t")
            ac = items[0]
            if ac not in all_ac2:
                all_ac2.append(ac)
            for j in items:
                #print j
                if j[0:3] == "IPR":
                    iprc2 += 1
                    if ac not in iacs2:
                        iacs2.append(ac)
                    #print j.split("; ")
                    for i in j.split("; "):
                        iprs2.append(i)
                        if i not in iprs22:
                            iprs22.append(i)
                elif j[0:3] == "GO:":
                    goc2 += 1
                    #print j.split("; ")
                    if ac not in gacs2:
                        gacs2.append(ac)
                    for i in j.split("; "):
                        gos2.append(i)
                        if i not in gos22:
                            gos22.append(i)
                
                

    

    ##############3   VENN DIAGGRAM   IPR1 vs. IPR2 #######################3

    print "TOTAL IPRs1 : ",len(iprs1),"\nTOTAL IPRs2 : " , len(iprs2)
    print "TOTAL GOs1 : ",len(gos1),"\nTOTAL GOs2 : " , len(gos2)

    print "UNIQUE IPRs1 : ",len(iprs11),"\nUNIQUE IPRs2 : " , len(iprs22)
    print "UNIQUE GOs1 : ",len(gos11),"\nUNIQUE GOs2 : " , len(gos22)


    ################################    VENN    ###################3
    one = set(iprs1) #uniq_iprs1)
    two = set(iprs2) # uniq_iprs2)

    # IPR IDs
    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", af, bf
    print "1 Union 2 = ", "\t", len(one.union(two))
    #a = most.common(one)
    
    print "\nDomains in common: ", af, bf
    print "1 Intersection 2 = ", "\t", len(one.intersection(two))
    #print one.intersection(two)

    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", af, bf
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - one.intersection(two))
    #**set(['d'])**

    #a = [1, 2, 3, 4, 5]
    #b = [9, 8, 7, 6, 5]
    #i for i, j in zip(a, b) if i == j]

    
    # GO TERMS
    three = set(gos1) #uniq_iprs1)
    four = set(gos2) # uniq_iprs2)

    #**set(['a', 'c', 'b'])**
    print "\nAll GOs ", af, bf
    print "2 Union 1 = ", "\t", len(four.union(three))
    #print four.union(three)
    
    print "\nGOs in common: ", af, bf
    print "2 Intersection 1 = ", "\t", len(four.intersection(three))
    #print one.intersection(two)

    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nGOs ", af, bf
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(four.union(three)  - four.intersection(three))

    ############ IPR SUM ##########
    
    print"\nTOTAL ACS 1 :",len(all_ac1),
    print"\tTOTAL ACS 2 :",len(all_ac2)
    
    c  = Counter(iprs1)
    print "\nTOTAL iprs1 found : ", len(iprs1) # 6597
    print "UNIQUE iprs1 : ", len(c), "\n", c.most_common(10) # 2088
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
    

    
    d  = Counter(iprs2)
    print "\nTOTAL IPRs2 found : ", len(iprs2) # 6597
    print "UNIQUE IPRs2 : ", len(d), "\n", d.most_common(10) # 3606
    #for i in dict(d.most_common()):
       # cf.write(str(i)+"\t"+str(d[i])+"\n")
        #print i, d[i]
    d_dict = dict(d.most_common())

    ### COVERAGE IPR 2
    #print iprc2, acs2
    try:
        cov = (len(iacs2)/float(len(all_ac2))) # most common(uncharacterized) / total found
        enr = (len(iprs22)/float(len(iprs2))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation IPR 2 :\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No IPRs 2 found!"
        #e  = Counter(ippr_names)
    #print "\nTOTAL IPRs NAMES found : ", len(ippr_names) # 6597
    #print "UNIQUE IPR NAMES : ", len(e), "\n", e.most_common(10) # 3606
    #for i in dict(e.most_common(10)):
     #   print i, e[i]


    ### GO SUM ###
    ### COVERAGE GO 1 
    
    f  = Counter(gos1)
    print "\nTOTAL GO IDS 1 found : ", len(gos1) # 6597
    print "UNIQUE GO IDs 1: ", len(f), "\n", f.most_common(10) # 1183
    #for i in dict(f.most_common()):
        #cf.write(str(i)+"\t"+str(f[i])+"\n")
        # print i, f[i]

    try:
        cov = (len(gacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enr = (len(gos11)/float(len(gos1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 1 :\n\tAc Coverage : ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No GOs 1 found!"
    

    g  = Counter(gos2)
    print "\nTOTAL GO 2 found : ", len(gos2)# 6597
    print "UNIQUE GO 2: ", len(g), "\n", g.most_common(10) # 1183
    g_dict = dict(g.most_common())

    
    ### COVERAGE GO 2
    try:
        cov = (len(gacs2)/float(len(all_ac2))) # most common(uncharacterized) / total found
        enr = (len(gos22)/float(len(gos2))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 2:\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No GOs 2 found!"

    print " Protein with IPRs 1: ", len(iacs1)
    print " Protein with IPRs 2: ", len(iacs2)

    print " Protein with GOs 1 : ", len(gacs1)
    print " Protein with GOs 2 : ", len(gacs2)


    #for i in g_dict:
        #cf.write(str(i)+"\t"+str(g[i])+"\n")
        #print i, g_dict[i]
        #oddsratio, pvalue = stats.fisher_exact([[1100,6848], [11860,75292]])

    ###############################    PVALUE + ODDSRATIO   ###################################

    cfile = str(af[:-4])+"_enrich_acgo.txt"
    cf = open(cfile, "w")
    dfile = str(af[:-4])+"_enrich_acipr.txt"
    df = open(dfile, "w")
    pvalues_ipr_sign = {}
    pvalues_go_sign = {}
    
    
    #print iprn_dict["IPR030470"]
    #iprn_dict["IPR030470"] = "null"
    df.write("IPR ID\tIPR Name\tcounts\ttotal_counts\tquery\tall\tp-value\todds_ratio\n") # lack id name. how to? fromgooripr_insert bmq_ipr2go
    results = []
    ipr_pv_dict = {}
    pvalues_ipr = []
    rank = 0
    pvalues = []
    for j in sorted(dict(c.most_common()), key=lambda x :dict(c.most_common())[x], reverse=True):
        #print j, f[j]
        if j in d_dict.keys():
            #print j, g_dict[j]
            n1 = c[j]
            n2 = d_dict[j]
            #n3 = len(c)- n1
            #n3 = len(iprs1) - n1  # n4 = len(d_dict) - n2
            
            #n4 = len(iprs2) - n2  # print j, n1, n2, n3, n4
            n3 = len(iacs1) - n1
            n4 = len(iacs2) - n2
            
            oddsratio, pvalue = stats.fisher_exact([[n1,n2], [n3,n4]],alternative="greater")
            #print oddsratio, pvalue
            #print '%f  ' % oddsratio
            #print '%f ' % pvalue
            

            pvalue = "{0:.10f}".format(pvalue)
            oddsratio = "{0:.10f}".format(oddsratio)

            #pvalues.append(pvalue)
            
                

            #print pvalue, oddsratio
            #oddsratio = remove_exponent(oddsratio)
            #pvalue = remove_exponent(pvalue)
            #print oddsratio, pvalue
            #if n1 and n2 and n3 and n4 > 0 else "null"

            #import numpy as np
            #import scipy.stats as stats

            #G = stats.poisson(30)
            #print(G.dist.expect(lambda x: (x+1), G.args, lb=0, ub=np.inf))
            # 31.0
            
            if j == "":
                pass
            if j not in iprn_dict.keys():
                iprn_dict[j] = "null"
            pvalues_ipr.append(pvalue)
            ipr_pv_dict[j]=pvalue
            
            #pvalue = lambda pvalue: float(pvalue.replace('.',','))
            #oddsratio = oddsratio = lambda oddsratio: float(oddsratio.replace('.',',')) 
            
            #print j, n1 , n2, n3, n4,"\tp-value = ", pvalue, "\tOdds Ratio = ", oddsratio
            #pvalue = float(pvalue)
            #oddsratio = float(oddsratio)
            
            
            #pvalue = abs(pvalue)
            #oddsratio = abs(oddsratio)
            #oddsratio = 1/float(oddsratio)
            
            results.append(str(pvalue)+"\t"+str(oddsratio))
            df.write(str(j)+"\t"+str(iprn_dict[j])+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)
                     +"\t"+str(pvalue)+"\t"+str(oddsratio)+"\n")

            if float(pvalue) <0.05:
                if j not in gon_dict.keys():
                    next
                pvalues_ipr_sign[j] = pvalue
                print str(j)+"\t"+str(iprn_dict[j])+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)+"\t"+str(pvalue)+"\t"+str(oddsratio)


    print "\nTOTAL Significant IPRs : ", len(pvalues_ipr_sign.keys()),"\n\n"
    
    cf.write("GO ID\tGO Term\tcounts\ttotal_counts\tquery\tall\tp-value\todds_ratio\n") # lack id name. how to? fromgooripr_insert bmq_ipr2go
    results = []
    go_pv_dict = {}
    pvalues_go = []
    for j in sorted(dict(f.most_common()), key=lambda x :dict(f.most_common())[x], reverse=True):
        #print j, f[j]
        if j in g_dict.keys():
            #print j, g_dict[j]
            n1 = f[j]
            n2 = g_dict[j]
            #n3 = len(gos1) - n1  # all gos
            #n4 = len(gos2) - n2 #
            n3 = len(gacs1) - n1  # all gos
            n4 = len(gacs2) - n2 # 
            
            #print j,f[j], g_dict[j], len(f), len(g_dict)
            #print n1, n2, n3, n4
            oddsratio, pvalue = stats.fisher_exact([[n1,n2], [n3,n4]],alternative="greater")
            #print float(oddsratio), float(pvalue)

            # "%.16f" % f if f >= 1e-16 else "0.0"
            pvalue = "{0:.10f}".format(pvalue)
            oddsratio = "{0:.10f}".format(oddsratio)
            #result_arr = fisherExact(np.array([[5, 0], [1, 4]]))
            #print result_arr
            
            # SCIPY pvalue, uses hypergeometric distribution (not efficient)
            #x = [[n1,n2],[n3,n4]]
            #result_pv = sp.stats.fisher_exact(x)
            # print result_pv[1], "pv"
            #pvalue2 = result_pv[1] # for i in result_pv:
            #odds2 = result_pv[0]
            
            #print float(oddsratio), float(pvalue)


            if j == "":
                pass
            if j not in gon_dict.keys():
                gon_dict[j] = "null"
                #break

            #pval_go = myformat(pvalue)
            pvalues_go.append(pvalue)
            go_pv_dict[j] = pvalue

            #pvalue = lambda pvalue: float(pvalue.replace('.',',')) # convertfloat = lambda x: float(str(x).replace(',','.')) 
            #oddsratio = lambda oddsratio: float(oddsratio.replace('.',',')) # convertfloat = lambda x: float(str(x).replace(',','.')) 
            
            #print j, n1 , n2, n3, n4,"\tp-value = ", pvalue, "\tOdds Ratio = ", oddsratio

            
            #pvalue = abs(pvalue)
            #oddsratio = abs(oddsratio)
            #oddsratio = 1/float(oddsratio)
            results.append(str(pvalue)+"\t"+str(oddsratio))
            cf.write(str(j)+"\t"+str(gon_dict[j])+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)+"\t"+str(pvalue)+"\t"+str(oddsratio)+"\n") # "\t"+str(pvalue2)+"\t"+str(odds2)+

            if float(pvalue) <0.05:
                if j not in gon_dict.keys():
                    print str(j)+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)+"\t"+str(pvalue)+"\t"+str(oddsratio)

                pvalues_go_sign[j] = pvalue
                print str(j)+"\t"+str(gon_dict[j])+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)+"\t"+str(pvalue)+"\t"+str(oddsratio)


    print "\nTOTAL Significant GOs : ", len(pvalues_go_sign.keys()), "\n"
   
    
        
    
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


    cf.close()
    df.close()


    # Plot 2y Scatter plot : ID vs p-value
    # Some points are missing while running the plot_2ylines riht away, so
    # we will split it in another function -> Scatter Enrich

    #print cfile, dfile
    #plot_2ylines2(cfile)
    #plot_2ylines2(dfile)
    #global theInputFiles
    theInputFiles.append(cf.name)
    theInputFiles.append(df.name) # NOT WORKING THE UPDATE
    #theFilenames1.append(cf.name)
    #global theFilenames1
    #theFilenames1.append(cf.name)
    #theFilenames1.append(df.name)
    
    #updatecombolist2(theInputFiles)
    

    # False Discovery Rate of the p-value 
        

    print cf.name, " saved!"
    print df.name, " saved!"
    #print ef, " saved!"
    #ef.close()
    #cff.close()
    #dff.close()
    kfile.close()
    #name1 = os.path.basename(cf)
    #name2 = os.path.basename(df)

    #global theInputFiles
    #theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)
    #Freq = 1000 # Set Frequency To 2500 Hertz
    #Dur = 100 # Set Duration To 1000 ms == 1 second
    #winsound.Beep(Freq,Dur)
    print "(complete!)"
    
#enrichment_analysis_domains("2pep_246_down_246_match_tax5671_bmq_ac_ipr_go.txt","output/tax5671_bmq_ac_ipr_go.txt","output/_bmq_ipr_go.txt")
#2pep_260_260_match_tax5671_bmq_ac_ipr_go.txt","output/tax5671_bmq_ac_ipr_go.txt","output/_bmq_ipr_go.txt")

#A4I5Q3_3990_bs_uniprot
# 246 vs all
#output/A4IBW8_entry_ipr_name_gene_go_path_local_exi.txt output/A4HU44_entry_ipr_name_gene_go_path_local_exi.txt

#enrichment_analysis_domains("output/2pep_246_down_filter_246_tax5671_bs_uniprot.txt","output/tax5671_bs_uniprot.txt", "output/_bmq_ipr_go.txt")



#enrichment_analysis_domains("2pep_246_down_246_match_tax5671_bmq_ac_ipr_go.txt","2peps4000_3998_match_tax5671_bmq_ac_ipr_go.txt","output/_bmq_ipr_go.txt")
#enrichment_analysis_domains("output/tax5671_bs_uniprot_filter_400_E9AG94.txt" , "output/tax5671_bs_uniprot_filter_5200_A4I2Z3.txt", "output/_bmq_ipr_go.txt")
"""
output/A4IBW8_entry_ipr_name_gene_go_path_local_exi.txt  New Selection 1 !
choice 1 =  output/A4IBW8_entry_ipr_name_gene_go_path_local_exi.txt
choice 2 =  output/A4I2I6_entry_ipr_name_gene_go_path_local_exi.txt
dict ipr2go =  output/_bmq_ipr_go.txt
TOTAL IPRs1 :  610 
TOTAL IPRs2 :  11239
TOTAL GOs1 :  252 
TOTAL GOs2 :  4971"""
#############################################################################################################################

def enrichment_analysis_accessions(af, bf, kf): # n = 1000

    """Using the most common dict (COLLECTIONS) file from last functions,
    we are able to build a default proteome set (bf) 
    and compare it with input proteome set (af) 
    combining de ipr2go EBI dictionary (kf)                                     """


    # Gets tax id and retrieves uni ids list from uniprot web
    # then enrichment analysis with list of accessions

    #from os.path import basename
    # now you can call it directly with basename
    #print basename(af)
    #str(af[:-4])+
    global theInputFiles
    
    print "choice 1 = ", af
    print "choice 2 = ", bf
    print "dict ipr2go = " , kf
    cfile = str(af[:-4])+"_enrich_tacsgo.txt"
    cf = open(cfile, "w")
    dfile = str(af[:-4])+"_enrich_tacsipr.txt"
    df = open(dfile, "w")
    
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
    num_lines1 = sum(1 for line in open(af))
    with open(af,"r") as f1:
        for i in f1.readlines():
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
                
                

    iacs2 = []
    #acs2 = []
    iprc2 = 0
    goc2 = 0
    iprs22 = []
    gos22 = []
    #num_lines2 = sum(1 for line in open(bf))
    with open(bf,"r") as f2:
        for i in f2.readlines():
            #print i
            i = i.strip()
            items = i.split("\t")
            ac = items[0]
            if ac not in all_ac2:
                all_ac2.append(ac)
            for j in items:
                #print j
                if j[0:3] == "IPR":
                    iprc2 += 1
                    if ac not in iacs2:
                        iacs2.append(ac)
                    #print j.split("; ")
                    for i in j.split("; "):
                        iprs2.append(i)
                        if i not in iprs22:
                            iprs22.append(i)
                elif j[0:3] == "GO:":
                    goc2 += 1
                    #print j.split("; ")
                    if ac not in gacs2:
                        gacs2.append(ac)
                    for i in j.split("; "):
                        gos2.append(i)
                        if i not in gos22:
                            gos22.append(i)
                
                

    

    ##############3   VENN DIAGGRAM   IPR1 vs. IPR2 #######################3

    print "TOTAL IPRs1 : ",len(iprs1),"\nTOTAL IPRs2 : " , len(iprs2)
    print "TOTAL GOs1 : ",len(gos1),"\nTOTAL GOs2 : " , len(gos2)

    print "UNIQUE IPRs1 : ",len(iprs11),"\nUNIQUE IPRs2 : " , len(iprs22)
    print "UNIQUE GOs1 : ",len(gos11),"\nUNIQUE GOs2 : " , len(gos22)


    ################################    VENN    ###################3
    one = set(iprs1) #uniq_iprs1)
    two = set(iprs2) # uniq_iprs2)

    # IPR IDs
    #**set(['a', 'c', 'b'])**
    print "\nAll Domains ", af, bf
    print "1 Union 2 = ", "\t", len(one.union(two))
    #a = most.common(one)
    
    print "\nDomains in common: ", af, bf
    print "1 Intersection 2 = ", "\t", len(one.intersection(two))
    #print one.intersection(two)

    #**set(['a', 'c', 'b', 'd'])**
    print "\nDomains ", af, bf
    print "1 Union 2 - 1 Intersection 2 = ","\t", len(one.union(two)  - one.intersection(two))
    #**set(['d'])**

    #a = [1, 2, 3, 4, 5]
    #b = [9, 8, 7, 6, 5]
    #i for i, j in zip(a, b) if i == j]

    
    # GO TERMS
    three = set(gos1) #uniq_iprs1)
    four = set(gos2) # uniq_iprs2)

    #**set(['a', 'c', 'b'])**
    print "\nAll GOs ", af, bf
    print "2 Union 1 = ", "\t", len(four.union(three))
    #print four.union(three)
    
    print "\nGOs in common: ", af, bf
    print "2 Intersection 1 = ", "\t", len(four.intersection(three))
    #print one.intersection(two)

    
    #**set(['a', 'c', 'b', 'd'])**
    print "\nGOs ", af, bf
    print "2 Union 1 - 2 Intersection 1 = ","\t", len(four.union(three)  - four.intersection(three))

    ############ IPR SUM ##########
    
    print"\nTOTAL ACS 1 :",len(all_ac1),
    print"\tTOTAL ACS 2 :",len(all_ac2)
    
    c  = Counter(iprs1)
    print "\nTOTAL iprs1 found : ", len(iprs1) # 6597
    print "UNIQUE iprs1 : ", len(c), "\n", c.most_common(10) # 2088
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
    

    
    d  = Counter(iprs2)
    print "\nTOTAL IPRs2 found : ", len(iprs2) # 6597
    print "UNIQUE IPRs2 : ", len(d), "\n", d.most_common(10) # 3606
    #for i in dict(d.most_common()):
       # cf.write(str(i)+"\t"+str(d[i])+"\n")
        #print i, d[i]
    d_dict = dict(d.most_common())

    ### COVERAGE IPR 2
    #print iprc2, acs2
    try:
        cov = (len(iacs2)/float(len(all_ac2))) # most common(uncharacterized) / total found
        enr = (len(iprs22)/float(len(iprs2))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation IPR 2 :\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No IPRs 2 found!"
        #e  = Counter(ippr_names)
    #print "\nTOTAL IPRs NAMES found : ", len(ippr_names) # 6597
    #print "UNIQUE IPR NAMES : ", len(e), "\n", e.most_common(10) # 3606
    #for i in dict(e.most_common(10)):
     #   print i, e[i]


    ### GO SUM ###
    ### COVERAGE GO 1 
    
    f  = Counter(gos1)
    print "\nTOTAL GO IDS 1 found : ", len(gos1) # 6597
    print "UNIQUE GO IDs 1: ", len(f), "\n", f.most_common(10) # 1183
    #for i in dict(f.most_common()):
        #cf.write(str(i)+"\t"+str(f[i])+"\n")
        # print i, f[i]

    try:
        cov = (len(gacs1)/float(len(all_ac1))) # most common(uncharacterized) / total found
        enr = (len(gos11)/float(len(gos1))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 1 :\n\tAc Coverage : ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No GOs 1 found!"
    

    g  = Counter(gos2)
    print "\nTOTAL GO 2 found : ", len(gos2)# 6597
    print "UNIQUE GO 2: ", len(g), "\n", g.most_common(10) # 1183
    g_dict = dict(g.most_common())

    
    ### COVERAGE GO 2
    try:
        cov = (len(gacs2)/float(len(all_ac2))) # most common(uncharacterized) / total found
        enr = (len(gos22)/float(len(gos2))) # unique / total found
        #'Correct answers: {:.2%}'.format(points/total) >'Correct answers: 88.64%'
        print "**** Estimation GO 2:\n\tAc Coverage: ", ("{0:.2%}".format(cov))#
        print "\tEnrichment", ("{0:.2%}".format(enr))
    except:
        print "No GOs 2 found!"

    print "1 Protein with IPRs: ", len(iacs1)
    print "2 Protein with IPRs: ", len(iacs2)

    print "1 Protein with GOs: ", len(gacs1)
    print "1 Protein with GOs: ", len(gacs2)


    #for i in g_dict:
        #cf.write(str(i)+"\t"+str(g[i])+"\n")
        #print i, g_dict[i]
        #oddsratio, pvalue = stats.fisher_exact([[1100,6848], [11860,75292]])

    ### PVALUE + ODDSRATIO
    
    #print iprn_dict["IPR030470"]
    #iprn_dict["IPR030470"] = "null"
    df.write("IPR ID\tIPR Name\tcounts\ttotal_counts\tquery\tall\tp-value\todds_ratio\n") # lack id name. how to? fromgooripr_insert bmq_ipr2go
    results = []
    ipr_pv_dict = {}
    pvalues_ipr = []
    rank = 0
    for j in sorted(dict(c.most_common()), key=lambda x :dict(c.most_common())[x], reverse=True):
        #print j, f[j]
        if j in d_dict.keys():
            #print j, g_dict[j]
            n1 = c[j]
            n2 = d_dict[j]
            #n3 = len(c)- n1
            #n3 = len(iprs1) - n1  # n4 = len(d_dict) - n2
            
            #n4 = len(iprs2) - n2  # print j, n1, n2, n3, n4
            n3 = len(iacs1) 
            n4 = len(iacs2) 
            
            oddsratio, pvalue = stats.fisher_exact([[n1,n2], [n3,n4]])
            #print oddsratio, pvalue
            #print '%f  ' % oddsratio
            #print '%f ' % pvalue

            pvalue = "{0:.10f}".format(pvalue)
            oddsratio = "{0:.10f}".format(oddsratio)
            #print pvalue, oddsratio
            #oddsratio = remove_exponent(oddsratio)
            #pvalue = remove_exponent(pvalue)
            #print oddsratio, pvalue
            #if n1 and n2 and n3 and n4 > 0 else "null"

            #import numpy as np
            #import scipy.stats as stats

            #G = stats.poisson(30)
            #print(G.dist.expect(lambda x: (x+1), G.args, lb=0, ub=np.inf))
            # 31.0
            
            if j == "":
                pass
            if j not in iprn_dict.keys():
                iprn_dict[j] = "null"
            pvalues_ipr.append(pvalue)
            ipr_pv_dict[j]=pvalue
            
            #pvalue = lambda pvalue: float(pvalue.replace('.',','))
            #oddsratio = oddsratio = lambda oddsratio: float(oddsratio.replace('.',',')) 
            
            #print j, n1 , n2, n3, n4,"\tp-value = ", pvalue, "\tOdds Ratio = ", oddsratio
            #pvalue = float(pvalue)
            #oddsratio = float(oddsratio)
            
            
            #pvalue = abs(pvalue)
            #oddsratio = abs(oddsratio)
            #oddsratio = 1/float(oddsratio)
            
            results.append(str(pvalue)+"\t"+str(oddsratio))
            df.write(str(j)+"\t"+str(iprn_dict[j])+"\t"+str(n1)+"\t"+str(n2)+"\t"+str(n3)+"\t"+str(n4)
                     +"\t"+str(pvalue)+"\t"+str(oddsratio)+"\n")
            

    cf.write("GO ID\tGO Term\tcounts\ttotal_counts\tquery\tall\tp-value\todds_ratio\n") # lack id name. how to? fromgooripr_insert bmq_ipr2go
    results = []
    go_pv_dict = {}
    pvalues_go = []
    for j in sorted(dict(f.most_common()), key=lambda x :dict(f.most_common())[x], reverse=True):
        #print j, f[j]
        if j in g_dict.keys():
            #print j, g_dict[j]
            n1 = f[j]
            n2 = g_dict[j]
            #n3 = len(gos1) - n1  # all gos
            #n4 = len(gos2) - n2 #
            n3 = len(gacs1)   # all gos
            n4 = len(gacs2)  # 
            
            #print j,f[j], g_dict[j], len(f), len(g_dict)
            #print n1, n2, n3, n4
            oddsratio, pvalue = stats.fisher_exact([[n1,n2], [n3,n4]])
            #print float(oddsratio), float(pvalue)

            # "%.16f" % f if f >= 1e-16 else "0.0"
            pvalue = "{0:.10f}".format(pvalue)
            oddsratio = "{0:.10f}".format(oddsratio)
            #result_arr = fisherExact(np.array([[5, 0], [1, 4]]))
            #print result_arr
            
            # SCIPY pvalue, uses hypergeometric distribution (not efficient)
            #x = [[n1,n2],[n3,n4]]
            #result_pv = sp.stats.fisher_exact(x)
            # print result_pv[1], "pv"
            #pvalue2 = result_pv[1] # for i in result_pv:
            #odds2 = result_pv[0]
            
            #print float(oddsratio), float(pvalue)


            if j == "":
                pass
            if j not in gon_dict.keys():
                gon_dict[j] = "null"
                #break

            #pval_go = myformat(pvalue)
            pvalues_go.append(pvalue)
            go_pv_dict[j] = pvalue

            #pvalue = lambda pvalue: float(pvalue.replace('.',',')) # convertfloat = lambda x: float(str(x).replace(',','.')) 
            #oddsratio = lambda oddsratio: float(oddsratio.replace('.',',')) # convertfloat = lambda x: float(str(x).replace(',','.')) 
            
            #print j, n1 , n2, n3, n4,"\tp-value = ", pvalue, "\tOdds Ratio = ", oddsratio

            
            #pvalue = abs(pvalue)
            #oddsratio = abs(oddsratio)
            #oddsratio = 1/float(oddsratio)
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


    cf.close()
    df.close()


    # Plot 2y Scatter plot : ID vs p-value
    # Some points are missing while running the plot_2ylines riht away, so
    # we will split it in another function -> Scatter Enrich

    #print cfile, dfile
    #plot_2ylines2(cfile)
    #plot_2ylines2(dfile)
    #global theInputFiles
    theInputFiles.append(cf.name)
    theInputFiles.append(df.name) # NOT WORKING THE UPDATE
    #theFilenames1.append(cf.name)
    #global theFilenames1
    #theFilenames1.append(cf.name)
    #theFilenames1.append(df.name)
    
    #updatecombolist2(theInputFiles)
    

    # False Discovery Rate of the p-value 
        

    print cf.name, " saved!"
    print df.name, " saved!"
    #print ef, " saved!"
    #ef.close()
    #cff.close()
    #dff.close()
    kfile.close()
    #name1 = os.path.basename(cf)
    #name2 = os.path.basename(df)

    #global theInputFiles
    #theInputFiles.append(cf.name)
    #theInputFiles.append(df.name)
    #Freq = 1000 # Set Frequency To 2500 Hertz
    #Dur = 100 # Set Duration To 1000 ms == 1 second
    #winsound.Beep(Freq,Dur)
    print "(complete!)"
    

######################################################################3
###########################################################################
#### P VALUE CORRECTION - FWER + FDR
#####################################################################
# test 1 : http://stats.stackexchange.com/questions/870/multiple-hypothesis-testing-correction-with-benjamini-hochberg-p-values-or-q-va
# 
#######################################################################

def FDR(x):
    """
    Assumes a list or numpy array x which contains p-values for multiple tests
    Copied from p.adjust function from R  
    """
    o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
    ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
    q = sum([1.0/i for i in xrange(1,len(x)+1)])
    l = [q*len(x)/i*x[j] for i,j in zip(reversed(xrange(1,len(x)+1)),o)]
    l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
    return l


#### AUXILIAR FUNCTIONS P_VALUE CORRECTION
    
def correct_pvalues_in_multiple_testing(pvalues, correction_type = "Benjamini-Hochberg"):                
    """                                                                                                   
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """
    from numpy import array, empty
    #print pvalues
    pvs = sorted(pvalues)
    pvalues = array(pvs)
    #print pvalues
    n = float(pvalues.shape[0]) # len
    new_pvalues = empty(n)
    new_values = {}
    old_values = {}
    
    ############################################## Pvalue Correction , Lesack 2011
    ##############################################
    if correction_type == "BonfCorr":                       # workes for small samples
        n = len(pvalues)
        for i,j in enumerate(pvs):
            pcorr = float(j) * n # Pcorr = P * n :: 0.01 * 300 = 3
            #pcorr = min(pcorr, 1)
            #pcorr = max(pcorr, 0)            
            new_values[i] = float(pcorr)
            
            old_values[i] = float(j)
            #print i, pcorr
            
    elif correction_type == "BonfHolmCorr":        # we reject if alpha / (n-rank) > pvalue    :: 0.05/295=0.0001                                           
        n = len(pvalues)
        for i,j in enumerate(pvalues):
            pcorr = float(n-i+2) * float(j) # Pcorr = p * (n-k+1) :: 0.01*(300-5+1)=2.95
            pcorr = min(pcorr, 1)
            pcorr = max(pcorr, 0) 
            new_values[i] = float(pcorr)
            old_values[i] = float(j)
            #print i, pcorr
    
    elif correction_type == "BenjHochCorr":                                                         
        n = len(pvalues)
        for i,j in enumerate(pvalues):
            pcorr = (float(j) * float(n))/(i+1)  # Pcorr = P*n / k :: 0.01*300/5 = 0.6:: 0.01*10/5=0.02

            #alpha = 0.05
            #pcorr = float(i+2/len(pvalues)) *float(alpha)
            new_values[i] = float(pcorr)
            old_values[i] = float(j)
            #print i, pcorr


    #############################################
    #############################################
            
    elif correction_type == "Benjamini-Hochberg":                                                         
        values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
        values.sort()
        #print values
        values.reverse()                                                                                  
        new_values = []
        for i, vals in enumerate(values):                                                                 
            rank = n - i
            pvalue, index = vals                                                                          
            new_values.append((n/rank) * float(pvalue))
            old_values[i] = vals
        for i in xrange(0, int(n)-1):  
            if new_values[i] < new_values[i+1]:                                                           
                new_values[i+1] = new_values[i]                                                         
        for i, vals in enumerate(values):
            pvalue, index = vals
            new_pvalues[index] = new_values[i]
            
        #print pvalues
        #print new_pvalues
        
    
    elif correction_type == "BjHc":                                                                   
        prev_bh_value = 0
        
        for i, p_value in enumerate(pvalues):
            bhc_value = float(p_value) * (len(pvalues) / (i+1))  # pv * m /(i)
            bh_value = float(p_value) * (len(pvalues) /(i+1)) #print p_value , bh_value
            bh_value = min(bh_value, 1)
            bh_value = max(bh_value, prev_bh_value)
            prev_bh_value = bh_value
            #print bh_value, bhc_value
            new_values[i] = bh_value
            old_values[i] = p_value

    ###############################      NEED REVIEW      ###################################
    elif correction_type == "BenjaminiHochbergCorr10":   # pcorr = p*len(p)/rank                                                              
        alpha = 0.10
        for i, p_value in enumerate(pvalues): # pv * (m/i) = adjusted
            #pcorr = float(p_value)* len(pvalues)/(i+1)

            bh_value = (float(p_value) *len(pvalues)) /(int(i+1))
            

            #print i, bh_value
            
            new_values[i] = bh_value
            old_values[i] = p_value
            
    elif correction_type == "BenjaminiHochbergCorr20":   # pcorr = p*len(p)/rank                                                              
        alpha = 0.20
        for i, p_value in enumerate(pvalues): # pv * (m/i) = adjusted
            #pcorr = float(p_value)* len(pvalues)/(i+1)

            bh_value = (float(p_value) *len(pvalues)) /(int(i+1))
            
            
            new_values[i] = bh_value
            old_values[i] = p_value
            
    elif correction_type == "BenjaminiHochbergCorr35":   # pcorr = p*len(p)/rank                                                              
        alpha = 0.35
        for i, p_value in enumerate(pvalues): # pv * (m/i) = adjusted
            #pcorr = float(p_value)* len(pvalues)/(i+1)

            bh_value = (float(p_value) *len(pvalues)) /(int(i+1))
            
            new_values[i] = bh_value
            old_values[i] = p_value

    ###################################################
    # TESTED for 3 alpha thresholds
    
    elif correction_type == "10BenjaminiHochberg":                                                                 
        alpha = 0.10
        for i, p_value in enumerate(pvalues): # pv * (m/i) = adjusted ### P < ( i/m ) Q
            
            bhh_value = (int(i+1)/ float(len(pvalues))) * float(alpha)    # 10/300 * 0.1 =0.003> 0.001 TRUE POSITIVE at FDR 10%
            bh1 = float(p_value) * (int(i+1)/ len(pvalues)) # 0.001 *10 /300 = 0.00003 NOT SIGNIFICANT
            pcorr = (float(p_value) * len(pvalues)) / (i+1) #print i, p_value
            #print int(i+1), len(pvalues), alpha
            
            
            new_values[i] = bhh_value
            old_values[i] = p_value

    elif correction_type == "20BenjaminiHochberg":                                                                   
        alpha = 0.20
        for i, p_value in enumerate(pvalues): # pv * (m/i) = adjusted  ## pcorr = p * n / rank
            bhh_value = (int(i+1)/ float(len(pvalues))) * float(alpha)
            
            bh1 = float(p_value) * (int(i+1)/ len(pvalues))
            #pcorr = (float(p_value) * len(pvalues)) / (float(i+1)
            
            new_values[i] = bhh_value
            old_values[i] = p_value

    elif correction_type == "35BenjaminiHochberg":                                                               
        alpha = 0.35
        for i, p_value in enumerate(pvalues): # pv * (m/i) = adjusted
            bhh_value = (int(i+1)/ float(len(pvalues))) * float(alpha)

            bh1 = float(p_value) * (int(i+1)/ len(pvalues))
            
            pcorr = (float(p_value) * len(pvalues)) / (i+1) 
            new_values[i] = bhh_value
            old_values[i] = p_value

    ########################################################################
    ########################################################################

            
    elif correction_type == "BonferroniHolm":  #    Bonferroni  Holm :  Step-down  :IN OUTPUT
        # m = len(pvals)
        # If the 1st p-value is greater than or equal to alpha/m, the procedure is stopped and no p-values are significant.
        #Otherwise, go on. 
        prev_bfhm_value = 0
        alpha = 0.10
        for i, p_value in enumerate(pvalues):
            #bfhm_value =  float(p_value) * (len(pvalues)- i + 2.0)
            #bfhm_value = (float(p_value) * len(pvalues))-i+1 # P <  m * p_value
            #print i, p_value, bfhm_value
            bfhm_value =  float(alpha) *(float(len(pvalues)- i))

            bfhm_value = min(bfhm_value, 1)
            bfhm_value = max(bfhm_value, prev_bfhm_value)
            
            new_values[i] = bfhm_value
            old_values[i] = p_value
            
    

    # TYPE I ERROR - Control False Positives - FWER     
    elif correction_type == "CalcBonf": # Pcorr = P * n                                                       
        
        for i, p_value in enumerate(pvalues): # pv * (m/i) = adjusted                      
            #print prev_bh_value
            pcorr = float(p_value) * len(pvalues)
            new_values[i] = pcorr
            old_values[i] = p_value
            
    elif correction_type == "CalcBonfHolm": # [ Pcorr = P * (n -k+1 ) ]
        
        for i, p_value in enumerate(pvalues): # pv * (m/i) = adjusted
            
            pcorr = float(p_value) * (len(pvalues)- i + 2.0)
            new_values[i] = pcorr
            old_values[i] = p_value

    # TYPE II ERROR - control False Negatives - FDR
    elif correction_type == "CalcBenjHoch": # [ Pcorr = P * n / k ]                                                    
        
       
        for i, p_value in enumerate(pvalues): # pv * (m/i) = adjusted
            
            pcorr = (float(p_value) * len(pvalues))/(i+1.0)
            new_values[i] = pcorr
            old_values[i] = p_value
    else:
        return
    
    return old_values, new_values

#print correct_pvalues_in_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1],correction_type = "Bonferroni")
#print correct_pvalues_in_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1],correction_type = "Bonferroni-Holm")
#print correct_pvalues_in_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1],correction_type = "Benjamini-Hochberg")
#correct_pvalues_in_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1],correction_type = "Bf")
#correct_pvalues_in_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1],correction_type = "BonferroniHolm")
#correct_pvalues_in_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1],correction_type = "BjHc")
#correct_pvalues_in_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1],correction_type = "BenjaminiHochberg")

#correct_pvalues_in_multiple_testing(pvalues, correction_type = "Benjamini-Holm")
#https://www.youtube.com/watch?v=cMFWFhTFohk

#############


#correct_pvalues_in_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1],correction_type = "Bonferroni-Holm")

###############################################################################3
##### CORRECT PVALES


def plot_2ylines_corr(df):

    # http://thomas-cokelaer.info/blog/2012/04/481/

    from collections import OrderedDict
    from collections import defaultdict
    import numpy as np

    from operator import itemgetter

    global theInputFiles
    #time_begin = datetime.datetime.fromtimestamp(time.time())

    

    # Open data to plot

    dfile = open(df, "r")
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
        elif "benjamini" == i:
            ind_bj = fields.index("benjamini")
            print i
        elif "bonferroni" == i:
            ind_bf = fields.index("bonferroni")
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
    #dict_ipr2pv = {}
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
        if len(items) == 8:
            #print items
            xlabel = items[0]
            #gos.append(items[0])
            lines.append(items)
            ipr = items[0]
            pv =items[ind_pv]
            #odds = i[ind_or]
            #bjmn = i[ind_bj]
            #bfrn = i[ind_bf]
            #pvals.append(items[-2])
            cts = items[2]
            #print cts, pv
            pvals6.append(pv)
            counts6.append(cts)
            #bjmn6.append(bjmn)
            #bfrn6.append(bfrn)
           
            #print cts
            #counts.append(items[2])
            #count2pv_dict[cts].append(pv) #= pv
            count2pv_dict[pv] = cts# append(pv)   USING DICT of P_VALUES
            joint = "\t".join(items)
            
            lines_id[ipr]= joint
            #print joint
            dict_ipr2pv[ipr] = str(pv)
            dict_ipr2cts[ipr] = str(cts)
            #dict_ipr2or[ipr] = str(odds)
            
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
        #bjmn.append(bjmn)
        #bfrn.append(bfrn)
           

    
    pvals2 = []
    counts2 = []
    
    lf = open(str(dfile.name[:-4])+"_corr_pvals.txt","w")
    header1 = "\t".join(fields[:-1])
    #lf.write(header1+"\trank\tfdr10\tfdr20\tfwer\n")
    lf.write(header1+"\trank\tBenjamini-Hockberg\tBonferroni-Holm\tBonferroni\n")

    #mf = open(str(dfile.name[:-4])+"_corr_or.txt","w")
    
    #mf.write(header1+"\trank\tfdr10\tfdr20\tfwer\n")

    
    pvals5 = []
    counts5 = []

    pvals4 = []
    counts4 = []
   
    

    
    order_by_pv = sorted(lines, key=itemgetter(ind_pv)) #lines.sort(key=lambda x: x[-2]) :error
    print len(order_by_pv), len(lines) #dict_ipr2pv
    
    pvals3 = []
    counts3 = []
    #pvs_sorted = OrderedDict(sorted(dict_ipr2pv.items(), key=lambda x: x[1]))

    pvs_sorted = sorted(dict_ipr2pv.iteritems(), key=lambda x:x[1]) # CUTS
    #print dict_ipr2pv
    for i, j in pvs_sorted: # ipr : pvalue : counts
        #print i, j , lines_id[i] # dict_ipr2cts[i]
        pvals3.append(j)
        counts3.append(dict_ipr2cts[i])
        #print i, dict_ipr2cts[i]

        
    # CORRECT VERSION
    
    bjmn2 = []
    bfrn2 = []
    bjmn22 = []
    
    #oldpvbf, newpvbf = correct_pvalues_in_multiple_testing(pvals6,correction_type = "Bf")
    #oldpvbfh, newpvbfh = correct_pvalues_in_multiple_testing(pvals6,correction_type = "Bonferroni-Holm")
    #oldpvbfh, newpvbfh = correct_pvalues_in_multiple_testing(pvals6,correction_type = "BonferroniHolm")
    #oldpvbjh, newpvbjh = correct_pvalues_in_multiple_testing(pvals6,correction_type = "BjHc") # out

    #oldpvbfc1, newpvbjhc1 = correct_pvalues_in_multiple_testing(pvals6,correction_type = "BonferroniCorr")
    #oldpvbfhc2, newpvbjhc2 = correct_pvalues_in_multiple_testing(pvals6,correction_type = "Bonferroni-HolmCorr")
    #oldpvbjhc3, newpvbjh3 = correct_pvalues_in_multiple_testing(pvals6,correction_type = "Benjamini-HochbergCorr") # out
    

    oldpvbjhc, newpvbjhc = correct_pvalues_in_multiple_testing(pvals6,correction_type = "BenjHochCorr") #"CalcBenjHoch")
    oldpvbfh, newpvbfh = correct_pvalues_in_multiple_testing(pvals6,correction_type = "BonfHolmCorr")# "CalcBonfHolm")
    oldpvbf, newpvbf = correct_pvalues_in_multiple_testing(pvals6,correction_type = "BonfCorr")#"CalcBonf")
    
    
    rank = 0
    for j ,i in enumerate(order_by_pv): # list inside list
        #print i
        if len(i) == 8:
            val = i[2]  #else "999"
            #print val
            new_line2 = "\t".join(i[:-1])
            
            
            #lf.write(new_line2+"\t"+str(j+1)+"\t"+str(newpvbjhc1[j])+"\t"+str(newpvbjhc2[j])+"\t"+str(newpvbfh[j])+"\n")
            lf.write(new_line2+"\t"+str(j+1)+"\t"+str(newpvbjhc[j])+"\t"+str(newpvbfh[j])+"\t"+str(newpvbf[j])+"\n")
            
            if int(val) < 1000:
                pvals2.append(i[ind_pv])
                counts2.append(i[2])
                bjmn2.append(newpvbjhc[j])
                bjmn22.append(newpvbfh[j])
                bfrn2.append(newpvbf[j])
                
                #print i[-2], i[2]
                
            else:
                print "hello i[2],i[ind_pv], val "
                print i[2],i[ind_pv], val
        else:
            return

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

    title = fields[0] +" vs. " +fields[ind_pv] + "\n" + (dfile.name)
    # create the general figure
    #fig1 = figure()    
    fig1 = plt.figure()
    fig1.canvas.set_window_title(str(lf.name[:-4]))

    fig1.suptitle(str(title))
    

    #t = np.arange(0.0, 1.0, 0.01)
    #s = np.sin(2*np.pi*t)
    #line, = ax.plot(t, s, color='blue', lw=2)
    
    

    # and the first axes using subplot populated with data
    ax1 = fig1.add_subplot(111)
    line1 = ax1.plot(pvals2, 'o-', label='p-value')
    #line3 = ax1.plot(bjmn2, 'o-', label='fdr .10')
    #line5 = ax1.plot(bjmn22,'o-', label='fdr .20')
    #line4 = ax1.plot(bfrn2, 'o-', label='fwer')

    line3 = ax1.plot(bjmn2, 'o-', label='fdr BHoch')
    line5 = ax1.plot(bjmn22,'o-', label='fwer BHolm')
    line4 = ax1.plot(bfrn2, 'o-', label='fwer Bf')
    
    #ax1.yaxis.tick_right()
    #ax1.yaxis.set_label_position("right")
    ylabel("P-value")
    legend(loc="upper left")
    #legend()

    # draw vertical line from (70,100) to (70, 250)
    ax1.plot([0, len(counts2)], [0.05, 0.05], 'g-', lw=2)

    #line0 = ax1.plot(t, bline, color='blue', lw=2)

    # now, the second axes that shares the x-axis with the ax1

    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
    line2 = ax2.plot(counts2, 'xr-', label='counts')
    
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ylabel("Counts")
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
    legend()
    #legend((line1, line2), ("p-value", "counts"))
    #host.xaxis_date()
    #fig1.show()
    #draw()
    show()
    
    #Freq = 1000 # Set Frequency To 2500 Hertz
    #Dur = 100 # Set Duration To 1000 ms == 1 second
    #winsound.Beep(Freq,Dur)
    #print "(complete!)"
        
    #time_end = datetime.datetime.fromtimestamp(time.time())
    #print("Time elapsed: ", str(time_end - time_begin))

#plot_2ylines_corr("20150924_enrch_greater/tax5671_bs_uniprot_filter_729_A4HWK0_enrich_acipr_order_pv.txt")

 
####################################################################################3
#####  ENRICHMENT ANALYSIS - IPR AND GO - CORRECTION
#####################################################################################3

import scipy as sp
import scipy.stats
import numpy as np

# 347515 L major
def plot_2ylines_fdr(df):

    # http://thomas-cokelaer.info/blog/2012/04/481/

    from collections import OrderedDict
    from collections import defaultdict
    import numpy as np

    from operator import itemgetter

    global theInputFiles
    #time_begin = datetime.datetime.fromtimestamp(time.time())

    

    # Open data to plot

    dfile = open(df, "r")
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
        elif "benjamini" == i:
            ind_bj = fields.index("benjamini")
            print i
        elif "bonferroni" == i:
            ind_bf = fields.index("bonferroni")
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
    #dict_ipr2pv = {}
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
        if len(items) == 8:
            #print items
            xlabel = items[0]
            #gos.append(items[0])
            lines.append(items)
            ipr = items[0]
            pv =items[ind_pv]
            #odds = i[ind_or]
            #bjmn = i[ind_bj]
            #bfrn = i[ind_bf]
            #pvals.append(items[-2])
            cts = items[2]
            #print cts, pv
            pvals6.append(pv)
            counts6.append(cts)
            #bjmn6.append(bjmn)
            #bfrn6.append(bfrn)
           
            #print cts
            #counts.append(items[2])
            #count2pv_dict[cts].append(pv) #= pv
            count2pv_dict[pv] = cts# append(pv)   USING DICT of P_VALUES
            joint = "\t".join(items)
            
            lines_id[ipr]= joint
            #print joint
            dict_ipr2pv[ipr] = str(pv)
            dict_ipr2cts[ipr] = str(cts)
            #dict_ipr2or[ipr] = str(odds)
            
        else:
            #print len(items), items
            return
        
            

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
        #bjmn.append(bjmn)
        #bfrn.append(bfrn)
           

    #pv_sorted = OrderedDict(sorted(dict_ipr2pv.itens(), key=lambda: x: x[1]))
    #it = iter(sorted(d.iteritems()))
    #for key in sorted(dict_ipr2pv):
    #    print dict_ipr2pv[key]

    
    
    pvals2 = []
    counts2 = []
    
    lf = open(str(dfile.name[:-4])+"_corr_fdrlevels.txt","w")
    header1 = "\t".join(fields[:-1])
    lf.write(header1+"\trank\tfdr35\tfdr20\tfdr10\tfwer\n")
    #lf.write(header1+"\trank\tBenjamini-Hockberg\tBonferroni-Holm\tBonferroni\n")

    #mf = open(str(dfile.name[:-4])+"_corr_or.txt","w")
    
    #mf.write(header1+"\trank\tfdr10\tfdr20\tfwer\n")

    
    pvals5 = []
    counts5 = []

    pvals4 = []
    counts4 = []
   

    
    order_by_pv = sorted(lines, key=itemgetter(ind_pv)) #lines.sort(key=lambda x: x[-2]) :error
    print len(order_by_pv), len(lines) #dict_ipr2pv
    
    pvals3 = []
    counts3 = []
    #pvs_sorted = OrderedDict(sorted(dict_ipr2pv.items(), key=lambda x: x[1]))

    pvs_sorted = sorted(dict_ipr2pv.iteritems(), key=lambda x:x[1]) # CUTS
    #print dict_ipr2pv
    for i, j in pvs_sorted: # ipr : pvalue : counts
        #print i, j , lines_id[i] # dict_ipr2cts[i]
        pvals3.append(j)
        counts3.append(dict_ipr2cts[i])
        #print i, dict_ipr2cts[i]

        
    # CORRECT VERSION
    
    bjmn35 = []
    bjmn20 = []
    bjmn10 = []
    bfrn = []
    
    #oldpvbf, newpvbf = correct_pvalues_in_multiple_testing(pvals6,correction_type = "Bf")
    #oldpvbfh, newpvbfh = correct_pvalues_in_multiple_testing(pvals6,correction_type = "Bonferroni-Holm")
    #oldpvbfh, newpvbfh = correct_pvalues_in_multiple_testing(pvals6,correction_type = "BonferroniHolm")
    #oldpvbjh, newpvbjh = correct_pvalues_in_multiple_testing(pvals6,correction_type = "BjHc") # out

    oldpvbf, newpvbf = correct_pvalues_in_multiple_testing(pvals6,correction_type = "BonferroniHolm") # out
    oldpvbjhc10, newpvbjhc10 = correct_pvalues_in_multiple_testing(pvals6,correction_type = "10BenjaminiHochberg")
    oldpvbjhc20, newpvbjhc20 = correct_pvalues_in_multiple_testing(pvals6,correction_type = "20BenjaminiHochberg")
    oldpvbjhc35, newpvbjhc35 = correct_pvalues_in_multiple_testing(pvals6,correction_type = "35BenjaminiHochberg")

    #oldpvbjhc, newpvbjhc = correct_pvalues_in_multiple_testing(pvals6,correction_type = "CalcBenjHoch")
    #oldpvbfh, newpvbfh = correct_pvalues_in_multiple_testing(pvals6,correction_type = "CalcBonfHolm")
    #oldpvbf, newpvbf = correct_pvalues_in_multiple_testing(pvals6,correction_type = "CalcBonf")
    
    
    
    rank = 0
    for j ,i in enumerate(order_by_pv): # list inside list
        #print i
        if len(i) == 8:
            val = i[2]  #else "999"
            #print val
            new_line2 = "\t".join(i[:-1])
            
            
            #lf.write(new_line2+"\t"+str(j+1)+"\t"+str(newpvbjhc1[j])+"\t"+str(newpvbjhc2[j])+"\t"+str(newpvbfh[j])+"\n")
            lf.write(new_line2+"\t"+str(j+1)+"\t"+str(newpvbjhc35[j])+"\t"+
                     str(newpvbjhc20[j])+"\t"+str(newpvbjhc10[j])+"\t"+str(newpvbf[j])+"\n")
            
            if int(val) < 1000:
                pvals2.append(i[ind_pv])
                counts2.append(i[2])
                bjmn35.append(newpvbjhc35[j])
                bjmn20.append(newpvbjhc20[j])
                bjmn10.append(newpvbjhc10[j])
                bfrn.append(newpvbf[j])
                
                
                
                #print i[-2], i[2]
                
            else:
                print "hello i[2],i[ind_pv], val "
                print i[2],i[ind_pv], val
        else:
            return

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

    title = fields[0] +" vs. " +fields[ind_pv] + "\n" + (dfile.name)
    # create the general figure
    #fig1 = figure()    
    fig1 = plt.figure()
    fig1.canvas.set_window_title(str(lf.name[:-4]))

    fig1.suptitle(str(title))
    

    #t = np.arange(0.0, 1.0, 0.01)
    #s = np.sin(2*np.pi*t)
    #line, = ax.plot(t, s, color='blue', lw=2)
    
    

    # and the first axes using subplot populated with data
    ax1 = fig1.add_subplot(111)
    line1 = ax1.plot(pvals2, 'o-', label='p-value')
    line3 = ax1.plot(bjmn35, 'o-', label='fdr .35')
    line6 = ax1.plot(bjmn20,'o-', label='fdr .20')
    line5 = ax1.plot(bjmn10,'o-', label='fdr .10')
    line4 = ax1.plot(bfrn, 'o-', label='fwer')

    #line3 = ax1.plot(bjmn2, 'o-', label='fdr BH')
    #line5 = ax1.plot(bjmn22,'o-', label='fwer BFH')
    #line4 = ax1.plot(bfrn2, 'o-', label='fwer BF')
    
    #ax1.yaxis.tick_right()
    #ax1.yaxis.set_label_position("right")
    ylabel("P-value")
    legend(loc="upper left")
    #legend()

    # draw vertical line from (70,100) to (70, 250)
    #ax1.plot([0, len(counts2)], [0.05, 0.05], 'g-', lw=2)

    #line0 = ax1.plot(t, bline, color='blue', lw=2)

    # now, the second axes that shares the x-axis with the ax1

    ax2 = fig1.add_subplot(111, sharex=ax1, frameon=False)
    line2 = ax2.plot(counts2, 'xr-', label='counts')
    
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ylabel("Counts")
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
    legend()
    #legend((line1, line2), ("p-value", "counts"))
    #host.xaxis_date()
    #fig1.show()
    #draw()
    show()
 


#plot_2ylines_corr("20150904_benjamini2/up4fold_filter_733_all_fold_filter_5214_tax5671_bs_uniprot_enrich_acipr_order_pv.txt")
#plot_2ylines_fdr("20150911_up4FC_1/all_up_4FC_733_filter_732_tax5671_bs_uniprot_enrich_acipr_order_pv.txt")
#plot_2ylines_corr("20150911_up4FC_1/all_up_4FC_733_filter_732_tax5671_bs_uniprot_enrich_acipr_order_pv.txt")

#plot_2ylines_fdr("20150924_enrch_greater/tax5671_bs_uniprot_filter_729_A4HWK0_enrich_acipr_order_pv.txt")

#############################################################################3

    

   
time_end = datetime.datetime.fromtimestamp(time.time())
print("Time elapsed: ", str(time_end - time_begin))

