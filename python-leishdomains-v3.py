# -*- coding: cp1252 -*-

global _authors
_authors ="""
###############################################################################
######                          PROTEOMICSTATS                          #######
###############################################################################

This module contains code ideas from
Computational Biochemistry Lessons 2012 by PhD Antonio E. Ferreira

<http://fc.ul.pt>

This file is part of Proteome Analysis Software

#
#  Copyright (c) 2014-2015 - FCUL - Antonio E. Ferreira <aeferreira@fc.ul.pt>
#
#  File author(s):

#      Jose A. Mendes , DQB-FCUL, <jmend3z@gmail.com>

    >tk_tax_id_v154.py                                      Visual interface
    >bioservices_functions_tk.py                            Biodata
        uniprot
        interpro
        quickgo
        kegg
        pdb
        ebi
    >run_biopython_tk.py                                    Biodata
        swissprot
    >summary_functions_tk.py                                Aux functions
    >summary_plot_tk.py                                     Graphics plots
    >summary_enrich_tk.py                                   Statistical analysis




#  With Collaborations :

#  Vinicius Cogo, DI-FCUL ,  <vvcogo@gmail.com>
    >run_tk.py
    >bi_biopython_functions_tk.py
        bioblast


#  Renato Fialho, DQB-FCUL ,
    >correct_ids_from_tab.py

#  Daniel Vilar, DQB-FCUL , <danvjmusic@gmail.com>
    >filter_domains_by_family.py


### License: GNU GPLv3 http://www.gnu.org/licenses/gpl.html
### General Public License

source code vs. object code ?

################################################################################
"""

import datetime, time
time_begin = datetime.datetime.fromtimestamp(time.time())

import string
print (string.punctuation)

import sys, os

import fnmatch
import shutil


print ("FILE PATH:", sys.argv[0])
print ("FILE NAME:", os.path.basename(sys.argv[0]))

sys.path.append("C:\Python276\Lib\site-packages")

## myroot = "C:\Users\\mendes\\Dropbox\\Thesis\\aeferreira\\leishdomains"



########################        IMPORT MODULE FUNCTIONS       ############################

print ("""\nSpecial dependencies are:
    Tkinter, 
    biopython; bioservices; xmltramp; wsdl; soappy;
    numpy; scipy; matplotlib;
    Venn3;Venn """)

#import bio as Bio
#from bioservices_functions_tk import *

#from run_bioblast_tk import * # from run_biopython_tk import *
#from bi_biopython_functions_tk import *



#from summary_functions_tk import *
#from auxiliary_functions_tk import *


#from summary_plot_tk import *
#from summary_enrich_tk import *







############################################################
###SOURCES COURSES
#http://www.greenteapress.com/thinkpython/html/book020.html
#http://www.python-course.eu/tkinter_entry_widgets.php

##############################################################

"""Program Interface for Proteomic Data Analysis"""
import tkinter
from tkinter import *
#import ttk

import winsound
from glob import *

#import tkMessageBox
#from tkMessageBox import *

#import tkFileDialog
#from tkFileDialog import *
#from fs import *

import sys
#from Bio import Entrez
#from numpy import *


global theInputFiles
theInputFiles = []

global choice
global choice2


########################################################################

#print "ALL MODULES REQUIRED FOR PROTEOMIC STATS:"


modulenames = set(sys.modules)&set(globals())
print (modulenames)

allmodules = [sys.modules[name] for name in modulenames]
#print allmodules



######################################################################
###        TKINTER - WIDJETS                                      ####
######################################################################
# https://www.python.org/psf/trademarks/

jp = Tk()

frontname = u"ProteomicStats\u2122"
jp.title(frontname)


######################################################################
##########               BUTTON  FUNCTIONS              ##############
######################################################################


def faznada():
    showinfo("About", "ProteomicStats \n Set 2015 \n  Jose' A. S. Mendes\n v 1.0 ")


def get_out():
    if askyesno('Verify', 'Are you sure?'):
        showinfo('Warning', 'Bye')
        quit()
    else:
        showinfo("Warning", "Exit cancelled.")


def mostra_erro():
    showerror("Error", "well, it doesn' t work yet :)")


def show_file():
    nome = askopenfilename()
    if len(nome) > 0:
        showinfo("OK", "Selected file:\n"+nome)
    else:
        showinfo("Info", "No file selected.")


def save_data():
    "This functions takes combined KEYS and makes a unique accession table separayed values"
    if askyesno('Verify', 'Save data?'):
        showinfo('Verify', 'Saving...')
        #SystemExit("Saving..."
    else:
        showinfo("Verify", "No data saved.")


def quit():
    global jp
    jp.destroy()

