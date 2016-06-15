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
print string.punctuation

import sys, os

import fnmatch

print "FILE PATH:", sys.argv[0]
print "FILE NAME:", os.path.basename(sys.argv[0])

sys.path.append("C:\Python276\Lib\site-packages")



import shutil
myroot = "C:\Users\\mendes\\Dropbox\\Thesis\\aeferreira\\leishdomains"


########################        IMPORT MODULE FUNCTIONS       ############################

print """\nSpecial dependencies are:
    Tkinter, 
    biopython; bioservices; xmltramp; wsdl; soappy;
    numpy; scipy; matplotlib;
    Venn3;Venn 
    
    
"""

from bioservices_functions_tk import *

from run_bioblast_tk import * # from run_biopython_tk import *
from bi_biopython_functions_tk import *



from summary_functions_tk import *
from auxiliary_functions_tk import *


from summary_plot_tk import *
from summary_enrich_tk import *







############################################################
###SOURCES COURSES
#http://www.greenteapress.com/thinkpython/html/book020.html
#http://www.python-course.eu/tkinter_entry_widgets.php

##############################################################

"""Program Interface for Proteomic Data Analysis"""

from Tkinter import *
import ttk

import winsound
from glob import *

from tkMessageBox import *

import tkFileDialog
#from tkFileDialog import *
#from fs import *

import sys
from Bio import Entrez
#from numpy import *


global theInputFiles
theInputFiles = []

global choice
global choice2


########################################################################

#print "ALL MODULES REQUIRED FOR PROTEOMIC STATS:"


modulenames = set(sys.modules)&set(globals())
print modulenames

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


################################################################
##################      OPEN ONE FILE FROM DIALOG      #########
################################################################

filename_input = ""


def open_file():

    #del theInputFiles[:]
    from glob import glob
    filename_input = tkFileDialog.askopenfilename(title="Open File List", filetypes=[("File TXT",".txt"),("All files",".*")])

    name = os.path.basename(filename_input)
    names = glob('*.txt')+glob('*/*.txt')

    if len(filename_input) > 0:
        my_text_box.delete(1.0 , END)

        name = os.path.basename(filename_input)

        for f in names:
            if name == f:
                #print "Dir Path: ", f

                name = f
                aFile = open(filename_input)
                for line in aFile.readlines():
                    #print line
                    my_text_box.insert( END, line )
                    name = os.path.basename(filename_input)
                    #print aFile.name
                    theInputFiles.append(name) # instead get total path to file



            elif name in f:
                #print "Dir Path: ", f
                name = f
                aFile = open(filename_input)
                #showinfo("OK", "Filename :\n"+filename_input)
                for line in aFile.readlines():
                    #print line
                    my_text_box.insert( END, line )
                    name = os.path.basename(filename_input)
                    #print aFile.name
                    theInputFiles.append(name) # instead get total path to file
        showinfo("Upload File", "Filename :\n"+filename_input)

    else:
        showinfo("Warning", "No file was selected.")

###################################################################################################################3


def save_file_name():
    from tkFileDialog import asksaveasfilename
    from tkFileDialog import showerror
    #filesavebox(msg=None, title=None, default=None)
    filename_input = tkFileDialog.asksaveasfile(defaultextension=".txt", title="Save File List", filetypes=[("File TXT",".txt"),("All files",".*")])
    if len(filename_input) > 0:
        showinfo("OK", "Filename :\n"+filename_input)

    else:
        showinfo("Warning", "No file was selected.")



#####################################################################################################################
##############################################        UPLOAD    LIST     ###########################################################
######################################################################################
from glob import glob

filename_input2 = ""

def upload_file_list():
    #del theFilenames1[:]

    from glob import glob

    global theFilenames1

    filename_input2 = tkFileDialog.askopenfilename(title="Open File List", filetypes=[("File TXT",".txt"),("All files",".*")])
    #name = os.path.basename(filename_input2)
    #names = glob('*.txt')+glob('*/*.txt')

    names1 = glob('*.txt')+glob('*/*.txt') #+glob('*.xls')
    names2 = glob('*.xls')+glob('*/*.xls')
    #names2 = glob('*.txt')+glob('*/*.txt')
    

    

    name = os.path.abspath(filename_input2) # string buffer
    namebase = os.path.basename(name)
    namedir = os.path.dirname(name)
    

    folder = namedir.split("\\")[-1]

    namewdir =folder + "\\" + namebase
    
    

    

    if len(filename_input2) > 0:
        aFile = open(filename_input2)
        #name = os.path.basename(filename_input2)
        #names = glob('*.txt')+glob('*/*.txt')
        for f in names1:
            if namebase == f:
                #print "Dir Path: ", f
                namebase = f
                
                my_text_box.delete(1.0 , END)
                for line in aFile.readlines():
                    #print line
                    my_text_box.insert( END, line)
                #showinfo("OK", "Filename :\n"+filename_input2)
                showname = namebase

            elif namewdir == f:
                #print "Dir Path: ", f
                namewdir = f
                
                my_text_box.delete(1.0 , END)
                for line in aFile.readlines():
                    #print line
                    my_text_box.insert( END, line)
                #showinfo("OK", "Filename :\n"+filename_input2)
                showname = namewdir
        for f in names2:
            if namebase == f:
                #print "Dir Path: ", f
                namebase = f
            
                #showinfo("OK", "Filename :\n"+filename_input2)
                showname = namebase

            elif namewdir == f:
                #print "Dir Path: ", f
                namewdir = f
                
                #showinfo("OK", "Filename :\n"+filename_input2)
                showname = namewdir
                
        showinfo("Upoad File", "Filename :\n"+filename_input2)
        theFilenames1.append(showname)

        #print aFile.name
        #w3 = Label(frm2, text = name)
        #w3.pack(side=TOP)

        #files = tuple(theFilenames1)
        #combo1(frm2,files)
        global choice
        #global choice2

        print  "FILE SELECTED : ", showname
        updatecombolist(theFilenames1)
        theInputFiles.append(showname)
        #updatecombolist2theFilenames1)
        #global box_value
        #box_value = StringVar()
        #box_value = nameewdasedas
        #print choice # _name_ac1.txt
        cbox.set(choice)  # manually update the var...
        #cbox2.set(choice2)
        #print box_value # PY_VAR25
        #return choice

    else:
        showinfo("Warning", "No file was selected.")


##############################################################################################################


filename_input3 = ""

def upload_file_list2():

    #del theFilenames1[:]

    filename_input3 = tkFileDialog.askopenfilename(title="Open File List", filetypes=[("File TXT",".txt"),("File EXCEL",".xls"),("All files",".*")])    

   
    names1 = glob('*.txt')+glob('*/*.txt') #+glob('*.xls')
    names2 = glob('*.xls')+glob('*/*.xls')
    #names2 = glob('*.txt')+glob('*/*.txt')
    global theFilenames1
    

    

    name = os.path.abspath(filename_input3) # string buffer
    namebase = os.path.basename(name)
    namedir = os.path.dirname(name)
    

    folder = namedir.split("\\")[-1]

    namewdir =folder + "\\" + namebase
    

    if len(filename_input3) > 0:
        
        for i in names1:
            aFile = open(filename_input3)
            if namebase == i:
                #if name == f:
                #print "Dir Path: ", f
                
                my_text_box.delete(1.0 , END)
                for line in aFile.readlines():
                    #print line
                    my_text_box.insert( END, line )
                #showinfo("OK", "Filename :\n"+filename_input3)
                showname = namebase
                


            elif namewdir == i:
                #if name == f:
                #print "Dir Path: ", f
                
                my_text_box.delete(1.0 , END)
                for line in aFile.readlines():
                    #print line
                    my_text_box.insert( END, line )
                #showinfo("OK", "Filename :\n"+filename_input3)
                showname = namewdir
        for i in names2:
            if namebase == i:
                #if name == f:
                #print "Dir Path: ", f
                
                #showinfo("OK", "Filename :\n"+filename_input3)
                showname = namebase
                


            elif namewdir == i:
                
                #showinfo("OK", "Filename :\n"+filename_input3)
                showname = namewdir
                

        
        showinfo("OK", "Filename :\n"+filename_input3)
        theFilenames1.append(showname)
        global choice2
        print  "FILE SELECTED :", showname
        #updatecombolist(theFilenames1)
        updatecombolist2(theFilenames1)
        theInputFiles.append(showname)
        aFile.close()

        #global box_value


        #box_value = StringVar()
        #box_value = nameewdasedas
        #print choice # _name_ac1.txt
        #cbox.set(choice)  # manually update the var...
        cbox2.set(choice2)
        #print box_value # PY_VAR25
        #return choice

    else:
        showinfo("Warning", "No file was selected.")



##############################################################################
##############              PRINT AUX FILES                ###################
##############################################################################

def print_readme(theFilename):
    f = open(theFilename, "r")
    for i in f.readlines():
        print i
    f.close()

def print_panel(theFilename3):
    f = open(theFilename3, "r")
    for i in f.readlines():
        print i
    f.close()


def print_authors():
    
    global _authors
    print _authors
    




theFilename2 = "taxIdentifiers_ncbi.txt"

def print_readme2(theFilename2):
    #window = Gui()
    f = open(theFilename2, "r")
    for i in f.readlines():
        print i
    f.close()
    #window.mainloop()



##############################################################################################3
##########################       Convert Tax Id from Name       ####################################3
##################################################################################


from Bio import Entrez
Entrez.email = "jmend3z@gmail.com"

def get_tax_id(species):
    """to get data from ncbi taxomomy, we need to have the taxid. we can
    get that by passing the species name to esearch, which will return
    the tax id"""
    species = species.replace(' ', "+").strip()
    search = Entrez.esearch(term = species, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    return record['IdList'][0]

def get_tax_data(taxid):
    """once we have the taxid, we can fetch the record"""
    search = Entrez.efetch(id = taxid, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)


def get_tax_id_w_name():
    if not Entrez.email:
        print "you must add your email address"
        sys.exit(2)
    #list = ['Leishmania infantum JPCM5','Helicobacter pylori 26695', 'Thermotoga maritima MSB8', 'Deinococcus radiodurans R1', 'Treponema pallidum subsp. pallidum str. Nichols', 'Aquifex aeolicus VF5', 'Archaeoglobus fulgidus DSM 4304']
    #i = iter(list)
    #item = i.next()
    #for item in list:
        #print item
    name = raw_input("Enter Specie Name: ")
    taxid = get_tax_id(name)
    data = get_tax_data(taxid)
    lineage = {d['Rank']:d['ScientificName'] for d in
                data[0]['LineageEx'] if d['Rank'] in ['phylum']}
    print lineage
    print "Tax ID is: %s" % taxid, "\n"

##################################################################

import webbrowser

def from_uni_id_get_interpro_browser(theUniIDs):

    print "Query Limit = 10 proteins "
    for i in theUniIDs[0:10]:
        webbrowser.open("http://www.ebi.ac.uk/interpro/protein/"+str(i))

# from_uni_id_get_interpro_browser("P51587"):

def from_uni_id_find_similar_proteins(theUniIDs):

    print "Query Limit = 10 proteins "
    for i in theUniIDs[0:10]:
        webbrowser.open("http://www.ebi.ac.uk/interpro/protein/"+str(i)+"/similar-proteins")

def from_kegg_paths_show_genes(theKeggPaths):

    print "Query Limit = 10 proteins "
    for i in theKeggPaths[0:10]:
        webbrowser.open("http://www.kegg.jp/entry/"+str(i))


def open_ncbi_blast():
    print "NCBI Blast in Browser . . . "
    webbrowser.open("http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch&PROG_DEF=blastn&BLAST_PROG_DEF=megaBlast&BLAST_SPEC=blast2seq")

def sequence_search_egg():
    webbrowser.open("http://eggnogdb.embl.de/#/app/seqscan")

#######################################################################################

##############################################################################
########################   TOPMENUS    +   FUNCTIONS   +   SESSIONS
##############################################################################


global fnam
fnam = ""
e = StringVar()


def okk(root5,e):
    
    print "Session Tag Name is : ", e.get()
    
    global fnam
    fnam = e.get()

    root5.destroy()


def pop_up_insert_name():

    #top = Toplevel(jp)
    root5=Tk()
    
    #new_window = Toplevel(root)
    root5.title("Saving Session . . . ")

    Label(root5, text="Save Session as :  ", height = 4, width = 5).pack()

    e = Entry(root5)
    e.pack(padx=8)


    b = Button(root5, text="OK", command=lambda: okk(root5,e))
    b.pack(pady=25,padx=50)
    

    root5.wait_window()
    #return fnam
    # when mainloop is on, the remaining funtions don t run


##################################################################3


def enter_name():


    #root3 = Tk()
    e1 = StringVar()
    new_window = Toplevel(jp)
    Label(new_window, text="Session Name : ").grid(sticky=W,row=0)
    e1=Entry(new_window,width=40).grid(row=0,column=1,sticky=W)

    #Label(new_window, text="Value 2").grid(pady=20,sticky=W,row=1)
    #e2=Entry(new_window,width=20).grid(row=1,column=1,pady=20,sticky=W)
    #nn = e1.get()
    okbut= Button(new_window, text="Enter",command= lambda: ok("OK")).grid(column=4,row=0,pady=30)
    #cancel = Button(new_window,text="Cancel",command=lambda: callback("CANCEL")).grid(column=1,row=4,pady=30)
    #print e1.get()
    #new_window.destroy()
    #return nn


    #root2 = Tk()
    # https://mail.python.org/pipermail/tkinter-discuss/2008-August/001618.html
    #global selected
    #txt=Listbox(root2)

#####################################################################################################
#####################################################################################################

def new_session():

    

    #fnam = ""
    # http://tkinter.unpythonic.net/wiki/tkFileDialog
    
    global myroot 

    new_path = str(datetime.datetime.now())
    #mypath = "".join(new_path)
    #print new_path.split(".")

    #import glob, os
    #os.chdir("/mydir")
    #for file in glob.glob("*.txt"):
        #print(file)

    folder = tkFileDialog.askdirectory(initialdir='.')

    print "FOLDER : ", folder
    os.chdir(folder)
    fnames = glob("*.txt") + glob("*.xls")
    os.chdir(myroot)
    
    for i, f in enumerate(fnames):
        print i, f
        namebase = os.path.basename(f)
        namedir = folder.split("/")[-1]
        namewdir = namedir + "\\" + namebase

        theInputFiles.append(namewdir)
        
   

 

    

    showinfo(" New Session ", "Folder imported :\n"+namewdir)
    #print "\n", af.name, "saved!"
    print "\nNew Session  . . .", folder, "\n", "complete!"




################################################################################################################





def save_session(theInputFiles):

    
    global myroot 

    new_path = str(datetime.datetime.now())
    #mypath = "".join(new_path)
    #print new_path.split(".")

    mpath = new_path.split(".")[0] 
    mypath = mpath.split(" ")
    #print mypath
    n0 = mypath[0]
    n00 =mypath[1]
    n1 = n0.split("-")
    n11 = "".join(n1)
    n2 = n00.split(":")
    n22 = "".join(n2)
    n3 = n11 # +n22

   
    ###  OPEN POP UP BOX

    print "###################    Saving Session    ########################"
    #name= raw_input("\Enter Tag: ")

    #namex = enter_name() #
    fname = pop_up_insert_name()
    #print fname # = Nome, because no return
    
    #print str(fnam)
    
    path = n3 + "_" + str(fnam) # "_".join(mypath)
    nopath = n3+"_"
    
    if fnam == "":
        showinfo("Save session", "No Session Saved!")
        return
    
    print "Folder : ", path
    
    if not os.path.isdir(path):
        os.makedirs(path) # NEW DIRECTORY



    pname =  str(path)+"\\SESSION_"+str(fnam)+".txt"
    af = open(pname, "w")
    print "Session File: ",pname


    for j,i in enumerate(theInputFiles):
        #print "\n",i
        #fname = i.split("/")[-1]
        fname2 = os.path.basename(i)
        #print fname2
        #fname2 = i.split("\\")[-1]
        #print fname2
        
        #print str(path)+fname2
        predest = str(path)+"\\"+fname2
        root = sys.argv[0]
        rootpaths = root.split("\\")[:-1]
        rootpath = "\\".join(rootpaths)
        #source  = rootpath+"\\"+i
        source = os.path.abspath(i)
        dest = myroot+"\\"+predest

        print j, "Source:  " , source,"\n", j, "Destiny: ",dest

        shutil.copyfile(source, dest)

        af.write(predest+"\n")
        # os.rename(source, dest)
        # shutil.move("path/to/current/file.foo", "path/to/new/destination/for/file.foo")

    showinfo(" Save Session ", "Complete :\n"+path)
    print "\n", af.name, "saved!"
    print "Session Saving . . .", path, "\n", "complete!"



def load_session():
    #name = raw_input("Enter File Name: ")

    print "\n##########    LOADING SESSION      ############"

    

    from glob import glob
    filename_input = tkFileDialog.askopenfilename(title="Open File Session", filetypes=[("File TXT",".txt"),("All files",".*")])

    name = os.path.basename(filename_input)
    names = glob('*.txt')+glob('*/*.txt')

    trig = name.split("_")

    


    if filename_input == "":
        showinfo("Load session", "No Session Loaded!")
        return

    if trig =="session" or "SESSION":
        #my_text_box.delete(1.0 , END)

        
        aFile = open(filename_input)
        i = 0
        
        for line in aFile.readlines():
            line = line.strip()
            print i, line
            theInputFiles.append(line)
            #print lineaf = open(afile, "r")
            #if len(i.split("."))==2:

            # for f in names:
            #if name == f:
            i+=1
            
        aFile.close()
        showinfo(" Load Session ", "Complete :\n"+name)
        print "\nSession Loading . . . ", name ,"\n","Complete!"
    else:
        showinfo(" Load Session ", "No Session loaded . . . ")


############################################################################
#######      VENN 3
###########################################################################

global choicev1
global choicev2
global choicev3
ven1 = ""
ven2 = ""
ven3 = ""
venn_e1 = StringVar()
venn_e2 = StringVar()
venn_e3 = StringVar()

def newselection_v1( event):
    value_of_combo = cboxv1.get()
    print(value_of_combo), " #  >>>>>    venn 3 File (1) !"
    global choicev1
    choicev1 = cboxv1.get()
    

def newselection_v2( event):
    value_of_combo = cboxv2.get()
    print(value_of_combo), " #  >>>>>   venn 3 File (2) !"
    global choicev2
    choicev2 = cboxv2.get()

def newselection_v3( event):
    value_of_combo = cboxv3.get()
    print(value_of_combo), " #  >>>>>    venn 3 File (3) !"
    global choicev3
    choicev3 = cboxv3.get()


def updatecombolist_v1(theFilenames1):

    for i in theInputFiles:
        if i not in theFilenames1:
            theFilenames1.append(i)

    files = tuple(theFilenames1)
    #files = tuple(theInputFiles)
    

    global choicev1

    box_valuev1 = StringVar()

    if len(theFilenames1) == 0:
        choicev1 = "null"
    else:
        choicev1 = theFilenames1[-1]
        
        box_valuev1 = StringVar()
        box_valuev1 = choicev1
        #print "Seleced : ", choice,"\n"
        cbox.set(choicev1)
        #print box_value # PY_VAR19, PY_VAR28

    cbox['values'] = files
    cbox.bind("<<ComboboxSelected>>", newselection_v1)
    cbox.current()


def updatecombolist_v2(theInputFiles):

    

    files = tuple(theInputFiles)

    global choicev2

    box_valuev2 = StringVar()

    if len(theFilenames1) == 0:
        choicev2 = "null"
    else:
        choicev2 = theInputFiles[-1]
        
        box_valuev2 = StringVar()
        box_valuev2 = choicev2
        #print "Seleced : ", choice,"\n"
        cbox.set(choicev2)
        #print box_value # PY_VAR19, PY_VAR28

    cbox['values'] = files
    cbox.bind("<<ComboboxSelected>>", newselection_v2)
    cbox.current()

def updatecombolist_v3(theInputFiles):

    files = tuple(theInputFiles)

    global choicev3

    box_valuev3 = StringVar()

    if len(theInputFiles) == 0:
        choice_v3 = "null"
    else:
        choice_v3 = theInputFiles[-1]
        
        box_valuev3 = StringVar()
        box_valuev3 = choice_v3
        #print "Seleced : ", choice,"\n"
        cbox.set(choice_v3)
        #print box_value # PY_VAR19, PY_VAR28

    cbox['values'] = files
    cbox.bind("<<ComboboxSelected>>", newselection_v3)
    cbox.current()





def ok_venn(root6,choicev1,choicev2,choicev3):
    
    print "Venn3 Diagram : %s %s %s " % (choicev1,choicev2, choicev3)
    
    
    get_domain_venn_diagram1(choicev1,choicev2,choicev3)


    root6.destroy()


def pop_up_venn3_files():

    #top = Toplevel(jp)
    root6=Tk()
    
    #new_window = Toplevel(root)
    root6.title(" Venn 3 Diagram  ")

    Label(root6, text="Select files :  ", height = 6, width = 20).pack()

    

    box_valuev1 = StringVar()
    cbox1 = ttk.Combobox(root6, width = 65, textvariable=box_valuev1, postcommand = lambda: updatecombolist_v1(theFilenames1))#,state='readonly')
    cbox1.pack(side=TOP)


    box_valuev2 = StringVar()
    cbox2 = ttk.Combobox(root6, width = 65, textvariable=box_valuev2, postcommand = lambda: updatecombolist_v2(theFilenames1))#,state='readonly')
    cbox2.pack(side=TOP)
    
    box_valuev3 = StringVar()
    cbox3 = ttk.Combobox(root6, width = 65, textvariable=box_valuev3, postcommand = lambda: updatecombolist_v3(theFilenames1))#,state='readonly')
    cbox3.pack(side=TOP)




    b = Button(root6, text="OK", command=lambda: ok_venn(root6,choicev1,choicev2,choicev3))
    b.pack(pady=25,padx=50,side=BOTTOM)
    

    root6.wait_window()
    #return fnam
    # when mainloop is on, the remaining funtions don t run
    

###############     CONTROL PANEL - IMAGE      ###########################33


def show_control_panel():
    
    novi = Toplevel()
    novi.title(" Control Panel ")
    canvas = Canvas(novi, width = 830, height = 250)
    canvas.pack(expand = YES, fill = BOTH)
    gif1 = PhotoImage(file = 'tk_control_panel_1.gif') #image not visual
    canvas.create_image(50, 10, image = gif1, anchor = NW)
    #assigned the gif1 to the canvas object
    canvas.gif1 = gif1


print_authors

###############################################################################
###############################################################################
######################        MENU                                    #########
###############################################################################


topmenu  = Menu(jp)
filemenu1 = Menu(topmenu)
filemenu2 = Menu(topmenu)
filemenu03 = Menu(topmenu)
filemenu3 = Menu(topmenu)
filemenu4 = Menu(topmenu)
filemenu5 = Menu(topmenu)
helpmenu = Menu(topmenu)

# topmenu is set in a main Window
jp.config(menu=topmenu)

# filemenu e helpmenu s�o submenus do topmenu, com diferentes os

# TOPMENU 1
topmenu.add_cascade(label="File", menu=filemenu1)
filemenu1.add_command(label="Open...", command=upload_file_list)
filemenu1.add_command(label="New Session", command = new_session)
filemenu1.add_separator()

filemenu1.add_command(label="Save Session", command=lambda: save_session(theInputFiles))
filemenu1.add_command(label="Load Session", command=load_session)

filemenu1.add_separator()

filemenu1.add_command(label="Exit", command=get_out)

#
# New  # Create # Import # Export # Workspace
#

# TOPMENU 2

topmenu.add_cascade(label="Search", menu=filemenu2)
#filemenu2.add_command(label="Taxonomic identifiers", command=faznada)

filemenu2.add_command(label="Get Tax ID from name", command=get_tax_id_w_name) #print_panel((theFilename3)))
filemenu2.add_command(label="Find protein name in organism", command= search_protein_name_in_org_bs_uniprot) #print_panel((theFilename3)))
filemenu2.add_command(label="Look for Pathways in organism", command= from_name_get_kegg_paths_genes) #print_panel((theFilename3)))


#Button(frame01,text="Get Tax ID from name",  bd=3, command=get_tax_id_w_name).pack(side=LEFT) ###### aux Kegg function from_uni_id_get_kegg:lambda : from_name_get_tax_id(species_list)
#Button(frame01,text="Find protein name in organism",  bd=3, command=search_protein_name_in_org_bs_uniprot).pack(side=LEFT) # aux UniProt function
#Button(frame01,text="Look for Pathways in organism",  bd=3 ,command=from_name_get_kegg_paths_genes).pack(side=LEFT) # aux Kegg fucntion

#filemenu2.add_separator()


#filemenu2.add_separator()
#
#filemenu2.add_command(label="Exit", command=get_out)

# TOPMENU 3

def goa_downloads():
    webbrowser.open("http://www.ebi.ac.uk/GOA/downloads")

def quickgo_datasets():
    webbrowser.open("http://www.ebi.ac.uk/QuickGO/Dataset.html")

def ebi_services():
    webbrowser.open("http://www.ebi.ac.uk/services")

def amigo_tools():
    webbrowser.open("http://amigo.geneontology.org/amigo/software_list")

def ebi_downloads():
    webbrowser.open("https://www.ebi.ac.uk/interpro/download.html")

def open_string():
    webbrowser.open("http://string-db.org/newstring_cgi/show_input_page.pl?UserId=MRjCrJ09y1ji&sessionId=M2qWtvnACPad")

def open_expasy():
    webbrowser.open("http://www.expasy.org/")

def open_pathcommons():
    webbrowser.open("http://www.pathwaycommons.org/about/")

def open_ipfam():
    webbrowser.open("http://www.ipfam.org/")


##### SEPARATOR ####

def tax_help():
    webbrowser.open('http://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi')


def uni_map():
    webbrowser.open('http://www.uniprot.org/mapping/')

def open_biomart_download():
    webbrowser.open('http://www.biomart.org/')


def uni_ac_lists():
    webbrowser.open('http://www.uniprot.org/uniprot/')

def open_interproscan_download():
    webbrowser.open('http://www.ebi.ac.uk/interpro/search/sequence-search')

def open_expasy_translate():
    webbrowser.open('http://web.expasy.org/translate/')

def uniprot_tax_db():
    webbrowser.open("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions")


def open_cytoscape():
    webbrowser.open("http://js.cytoscape.org/")

def open_kegg_map_color():
    webbrowser.open("http://www.genome.jp/kegg/tool/map_pathway2.html")

def open_kegg_map_construct():
    webbrowser.open("http://www.genome.jp/kegg/tool/map_module.html")




#topmenu.add_cascade(label="Session", menu=filemenu03)
#filemenu03.add_command(label="Save Session", command=lambda: save_session(theInputFiles))
#filemenu03.add_command(label="Load Session", command=load_session)

#bt20=Button(frame033,text=" Save Session ", command=lambda: save_session(theInputFiles)).pack(padx=5, side=LEFT, anchor=N, expand=YES)# retrieve_input
#bt21=Button(frame033,text=" Load Session ", command= load_session).pack(padx=5, side=LEFT, anchor=N, expand=YES)# retrieve_input


topmenu.add_cascade(label="Links", menu=filemenu3)


# TOOLS
filemenu3.add_command(label="NCBI Tax ID Helper", command=tax_help)
filemenu3.add_command(label="UniProt Mapping", command=uni_map)
filemenu3.add_command(label="Expasy Translate DNA/RNA", command=open_expasy_translate)

filemenu3.add_command(label="BioMart Tools", command=open_biomart_download)

filemenu3.add_command(label="InterProScan", command=open_interproscan_download)

filemenu3.add_command(label="AmiGO Tools", command=amigo_tools)


filemenu3.add_separator()

filemenu3.add_command(label="EBI Services", command=ebi_services)
filemenu3.add_command(label="EBI Downloads", command=amigo_tools)

filemenu3.add_command(label="UniProt ID lists", command=uni_ac_lists)
filemenu3.add_command(label="UniProt Taxonomic Divisions ", command=uniprot_tax_db)


filemenu3.add_command(label="UniProt GOA Downloads", command=goa_downloads)

filemenu3.add_command(label="QuickGO Datasets", command=quickgo_datasets)

filemenu3.add_separator()

#LINKS

filemenu3.add_command(label="String 1.0 ", command = open_string)
filemenu3.add_command(label="ExPASy ", command = open_expasy)
filemenu3.add_command(label="Pathway Commons", command = open_pathcommons)
filemenu3.add_command(label="iPfam", command = open_ipfam)

filemenu3.add_separator()


filemenu3.add_command(label="Cytoscape", command = open_cytoscape)

filemenu3.add_command(label="Kegg Map Color", command = open_kegg_map_color)

filemenu3.add_command(label="Kegg Construct Module", comman= open_kegg_map_construct)





#Button(frame1, text='NCBI Tax ID Helper', fg= "blue", command=tax_help, cursor="hand2").pack(side=LEFT) # opens WEB

# Button(frame1,text="UniProt Mapping", fg="blue", command=uni_map, cursor="hand2").pack(side=LEFT) # opens WEB
# Button(frame1, text="UniProt ID lists", fg="blue", command=uni_ac_lists ,cursor="hand2").pack(side=LEFT)

# Button(frame1, text="BioMart Tools", fg="blue", command=open_biomart_download ,cursor="hand2").pack(side=LEFT)
# Button(frame1, text="InterProScan", fg="blue", command=open_interproscan_download ,cursor="hand2").pack(side=LEFT)


##########
# TOPMENU 5
# Ipython Examples
def ipynb_uniprot():
    webbrowser.open("http://nbviewer.ipython.org/url/pythonhosted.org//bioservices/_downloads/UniProt.ipynb")

topmenu.add_cascade(label="IPython notebooks", menu=filemenu4)
filemenu4.add_command(label="Uniprot example 1 ", command=ipynb_uniprot)

#filemenu.add_command(label="QuickGO Datasets", command=quickgo_datasets)
#filemenu.add_command(label="EBI Services", command=ebi_services)
#filemenu.add_command(label="AmiGO Tools", command=amigo_tools)

#filemenu.add_command(label="QuickGO Datasets", command=quickgo_datasets)


#filemenu.add_separator()

###################
# TOPMENU 6
# DOCUMENTATIONS : BioServices + BioPython

def bs_uniprot_doc():
    webbrowser.open("http://pythonhosted.org/bioservices/_modules/bioservices/uniprot.html")

def biopython_doc():
    webbrowser.open("http://biopython.org/wiki/Main_Page")

def fisher_link():
    webbrowser.open("http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.fisher_exact.html")

def correction_link():
    webbrowser.open("http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python")

   
def goatools_link():
    webbrowser.open("https://pypi.python.org/pypi/goatools")
##

def tk_tutorial():
    webbrowser.open("http://www.tkdocs.com/tutorial/")

def tkinter_wiki():
    webbrowser.open("https://wiki.python.org/moin/TkInter")

def lib_collections():
    webbrowser.open("https://docs.python.org/2/library/collections.html")
##


topmenu.add_cascade(label="Documentation ", menu=filemenu5)

filemenu5.add_command(label="BioServices UniProt Code", command=bs_uniprot_doc)

filemenu5.add_command(label="BioPython Wiki", command=biopython_doc)


filemenu5.add_separator()
filemenu5.add_command(label="Tkinter Tutorial", command=tk_tutorial)
filemenu5.add_command(label="Tkinter Wiki", command=tkinter_wiki)
filemenu5.add_command(label="Collections Lib", command=lib_collections)

filemenu5.add_separator()

filemenu5.add_command(label="Fisher Exact test Lib", command=fisher_link)

filemenu5.add_command(label="Multiple Tests Lib", command=correction_link)


filemenu5.add_separator()
filemenu5.add_command(label="GOA TOOLS Lib ", command=fisher_link)



##################
# TOPMENU 10
# Preferences # Workspace manager # Updates




theFilename3 = "quick_tutorial.txt"
theFilename4 = "enrich_tutorial.txt"

topmenu.add_cascade(label="Help", menu=helpmenu)
helpmenu.add_command(label="About...", command=faznada)

helpmenu.add_command(label="Tax Id list", command=lambda: print_readme2(theFilename2))

helpmenu.add_command(label="Quick Guide", command=lambda: open_list_box1 (theFilename3)) #print_panel((theFilename3)))
helpmenu.add_command(label="Enrichment Guide", command= lambda: open_list_box1(theFilename4)) #print_panel((theFilename3)))

helpmenu.add_command(label="Control Panel", command= show_control_panel) #print_panel((theFilename3)))

helpmenu.add_command(label="Authors", command = print_authors)




###############################################################################
##############################      FRAMES and BUTTONS     ###########################
#################################################################################

import webbrowser

# USEFUL TOOL IN UNIPROT
# http://www.ebi.ac.uk/uniprot/search/tools.html

# CREATE LABEL BUTTONS

#frontname = u"ProteomicStats\u2122"
lb = Label(jp)
lb.config(text = " ProteomicStats " , bg='black', fg='blue', font=('times', 20, 'bold'), height=2, width=30) # fam�lia, tamanho, estilo
lb.pack(expand=YES, fill=BOTH)
#l.grid(row=0, column = 0)


#####   INITIAL FRAME  #####

#frame1 = Frame(jp)# bg = "orange"
#frame1.pack(expand =YES, side=TOP)


###############

# LEISH BIOCYC
# http://www.biocyc.org/LEISH/organism-summary?object=LEISH



# USE EBI INTERPRO TO GET PROTEIN WITH GO IDS : Statistics
# http://www.ebi.ac.uk/QuickGO/GProtein?ac=P06727

# USE INTERPRO TO GET SEQ INFORMATION
# http://www.ebi.ac.uk/interpro/#

# ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/external2go/uniprotkb_sl2go
# ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/external2go/unipathway2go
# http://www.geneontology.org/external2go/ec2go

# http://www.genome.jp/dbget-bin/www_bget?hsa:7529

#Button(frame1, text="Tax ID help list", command= lambda: print_readme(theFilename2)).pack(side=LEFT)

### OPENS LINKS




# Button(frame1, text='NCBI Tax ID Helper', fg= "blue", command=tax_help, cursor="hand2").pack(side=LEFT) # opens WEB

# Button(frame1,text="UniProt Mapping", fg="blue", command=uni_map, cursor="hand2").pack(side=LEFT) # opens WEB
# Button(frame1, text="UniProt ID lists", fg="blue", command=uni_ac_lists ,cursor="hand2").pack(side=LEFT)

# Button(frame1, text="BioMart Tools", fg="blue", command=open_biomart_download ,cursor="hand2").pack(side=LEFT)
# Button(frame1, text="InterProScan", fg="blue", command=open_interproscan_download ,cursor="hand2").pack(side=LEFT)




frame01 = Frame(jp, bg = "white")
frame01.pack(expand =YES, side=TOP)

##Button(frame01,text="Match sequence",  bd=3, command=from_seq_get_bs_uniprot).pack(side=LEFT) # aux UniProt fucntion

# from_uni_ids_get_kegg_path
#Button(frame01,text="Look for Pathways",  bd=3 ,command=from_name_get_kegg_paths_genes).pack(side=LEFT) # aux Kegg fucntion

# http://www.ebi.ac.uk/interpro/search?q=kinase

#####   FRAME  #####

frame02 = Frame(jp, width=300,height=200)
frame02.pack(side=TOP,  expand = YES )


v1 = StringVar()
v2 = StringVar()
v3 = StringVar()
v4 = StringVar()
v5 = StringVar()
v6 = StringVar()

v1.set("5671")
v2.set("A4HUD2")#A4HUD2	A4HUD2_LEIIN	unreviewed	Uncharacterized protein	LINJ_10_0240		Predicted
v3.set("IPR002085") #'A1Y2C6 -> IPR002085 -> Alcohol dehydrogenase superfamily, zinc-type > GO:zinc ion binding > GO:0008270 > function\n
v4.set("GO:0008270")
v5.set("LINJ_10_0240")
#v6.set("7545")
v6.set("1FDH")



###########################################################################
###################################        CALLBACKS         ##############
def cb1():
    print "\tTax ID : %s" %E01.get()
    my_text_box.delete(1.0, END)
    my_text_box.insert( END,  E01.get()+"\n")

def cb2():
    print "\tUniProt ID : %s" %E02.get()
    my_text_box.delete(1.0, END)
    my_text_box.insert( END,  E02.get()+"\n")

def cb3():
    print "\tInterPro ID : %s" %E03.get()
    my_text_box.delete(1.0, END)
    my_text_box.insert( END,  E03.get()+"\n")

def cb4():
    print "\tGO ID : %s" %E04.get()
    my_text_box.delete(1.0, END)
    my_text_box.insert( END,  E04.get()+"\n")

def cb5():
    print "\tGene ID : %s" %E05.get()
    my_text_box.delete(1.0, END)
    my_text_box.insert( END,  E05.get()+"\n")

def cb6():
    print "\tPDB ID : %s" %E06.get()
    my_text_box.delete(1.0, END)
    my_text_box.insert( END,  E06.get()+"\n")


###########################    RADIOBUTTONS   ###########################3
# are associated to the same v, but different value=*

v = StringVar()
v.set("1") # inicializar a vari�vel com um val�r usando set()
           # sen�o, todos os bot�es ficam marcados.



b1 = Radiobutton(frame02, text="Enter Taxonomic ID",  command=cb1, bd=3,  variable=v, value="1")
b1.grid(row=1, column=0, sticky="W")
b2 = Radiobutton(frame02, text="Enter UniProt ID", command=cb2, bd=3,  variable=v, value="2")
b2.grid(row=1, column=1, sticky="W")
b3 = Radiobutton(frame02, text="Enter InterPro ID",  command=cb3, bd=3, variable=v, value="3")
b3.grid(row=1, column=2, sticky="W")
b4 = Radiobutton(frame02, text="Enter GO ID",  command=cb4, bd=3, variable=v, value="4")
b4.grid(row=1, column=3, sticky="W")
b5 = Radiobutton(frame02, text="Enter Gene/Kegg ID",  command=cb5, bd=3, variable=v, value="5")
b5.grid(row=1, column=4, sticky="W")
b6 = Radiobutton(frame02, text="Enter PDB ID",  command=cb6, bd=3, variable=v, value="6")
b6.grid(row=1, column=5, sticky="W")



E01 = Entry(frame02, textvariable=v1, bd=2)
E01.grid(row=2, column=0, sticky="W")
E01.focus()                                        # save a click
E01.bind('<Return>', (lambda event: cb1()))      # on enter key

E02 = Entry(frame02, textvariable=v2, bd=2)
E02.grid(row=2, column=1, sticky="W")

E03 = Entry(frame02, textvariable=v3, bd=2)
E03.grid(row=2, column=2, sticky="W")

E04 = Entry(frame02, textvariable=v4, bd=2)
E04.grid(row=2, column=3, sticky="W")

E05 = Entry(frame02, textvariable=v5, bd=2)
E05.grid(row=2, column=4, sticky="W")

E06 = Entry(frame02, textvariable=v6, bd=2)
E06.grid(row=2, column=5, sticky="W")

########################################################################3

################33
#frame022 = Frame(jp)
#frame022.pack()
#bt1 = Button(frame022, text='Upload File ... ', command=open_file)
#bt1.grid( row= 0, column=3, sticky="W")
#bt1.pack(side=LEFT)
#3bt2 = Label(frame022, text=filename_input)
#bt2.pack(side=RIGHT)
#bt1.grid( row=0, column=4, sticky="W")

########   FRAME   #############################33

frame033= Frame(jp)
frame033.pack()
bt10 =Button(frame033, text=" Upload List ", command=upload_file_list).pack(padx=5,side=LEFT, anchor=N, expand=YES)# open_file

#lb10 =Label(frame033, text=name).pack(side = RIGHT, fill=X)


frame03 = Frame(jp, height=30, width=250)
frame03.pack()


scrollbar = Scrollbar(frame03)
scrollbar.pack( side = LEFT, fill=Y )


my_text_box = Text(frame03, height=5, width=80,  wrap=WORD, yscrollcommand=scrollbar.set) # CHANGE TEXT BOX SIZE HERE

#my_text_box.grid(row=0, column=0, sticky=N+S+E+W)
my_text_box.pack(side = LEFT, fill = BOTH)

my_text_box.mark_set(INSERT, 1.0) #("insert", "%d.%d" % (line + 1, column + 1)


text1 = ["Insert","ID", "list", "here"]

for line in text1:
    my_text_box.insert( END, line +"\n")


def save_my_text():
    name = raw_input("Enter File Name: ")
    aFile = open("output/"+str(name)+".txt", "w")
    text1 =my_text_box.get('1.0',END).splitlines()
    for line in text1:
        aFile.write(line+"\n")
    aFile.close()
    theInputFiles.append(aFile.name)
    print aFile.name, "saved!"


def file_save_from_box():

    # http://stackoverflow.com/questions/19476232/save-file-dialog-in-tkinter

    mypath = "C:\Users\mendes\Dropbox\Thesis\aeferreira\leishdomains"
    mydir = "leishdomains"
    

    

    text1 = my_text_box.get('1.0',END)#.splitlines()
    text2 = []
    names0 = glob('*.txt') 
    names1 = glob('*/*.txt')
    #print names1
        
    f = tkFileDialog.asksaveasfile(title="Save File as ", mode='w', defaultextension=".txt")

    
    if f :
        for i in text1:
            f.write(i)
        f.close()
        #fullname = f.name # full name
        #namedir = fullname.split("\\")[-3:]
        #namedir2 = "\\".join(namedir)
        #name = os.path.basename(fullname)

        #name = nam.split(nam)[-1]
        #if name in names:
        
        namepath = os.path.abspath(f.name) # string buffer
        #namebase = os.path.basename(namepath)
        #namedir = os.path.dirname(namebase)
        #print namebase, namedir

        namebase = os.path.basename(f.name)
        namedir = os.path.dirname(f.name)
        #print namebase, namedir


        #name = os.path.splitext(f)[0] # ('file', '.ext')
        #infilename = namedir2 #+"\\"+ name
        folder = namedir.split("/")[-1]
        base = namebase.split("/")[-1]
        namewdir = folder + "\\\\" + base
        
        #print  namewdir #20150902_significantscat\\edfas.txt, still doesn t find
        showname = base

        #if base in names0:
            #showname = base
        #elif namewdir in names1:
            #showname = namewdir

        if folder == mydir:
            showname = base
        elif folder != mydir:
            showname = namewdir

        theInputFiles.append(showname)
        print "\n",showname, "saved!"
        
            
        
        showinfo(" File Saved ", "Filename :\n"+namepath)

    else:
        showinfo("Warning", "No file was saved.")# asksaveasfile return `None` if dialog closed with "cancel".
        #text2save = str(text1.get(1.0, END)) # starts from `1.0`, not `0.0`

#file_save()


bt19=Button(frame033,text=" Save List ", command=file_save_from_box).pack(padx=5, side=LEFT, anchor=N, expand=YES)# retrieve_input
#bt19=Button(frame033,text=" Save List ", command=save_file_name).pack(padx=5, side=LEFT, anchor=N, expand=YES)# retrieve_input


#A4IDS4
#########   NOTEPAD  ########

      

def opennotepadFile(choice):
    fileName = choice
    os.system("notepad.exe " + fileName)

def openexcelfile(choice):

    # sys.path[0]
    # change dir
    # os.chdir(sys.path[0])
    # or give full path for the file you want to open
    # os.system('start excel.exe "%s\\file.xls"' % (sys.path[0], ))
    print sys.path[0]
    mypath = os.path.abspath(choice)

    os.system('start excel.exe %s' % (mypath))



def open_file_edit(choice):

    ext = choice[-4:]
    if ext == ".txt":
        opennotepadFile(choice)
    elif ext == ".xls":
        openexcelfile(choice)
    else:
        return

#choice = "test34.xls"
#open_file_edit(choice)

#######################################################


#bt21=Button(frame033,text="Edit in Notepad ", command= lambda: opennotepadFile(choice)).pack(padx=5, side=LEFT, anchor=N, expand=YES)# retrieve_input



scrollbar.config( command = my_text_box.yview )
scrollbar.pack(side=LEFT)


ac_IDS = []
def retrieve_input():
    text1 =my_text_box.get('1.0',END).splitlines() # retrieve as list
    # input = self.myText_Box.get("1.0",'end-1c')
    print text1
    for i in text1:
        print i
        ac_IDS.append(i)
    print ac_IDS

####################################################################################
import io
def count_input():
    text1 = io.StringIO(my_text_box.get('1.0', END))
    IDs = 0
    for line in text1:
        line = line.split("\t")[0]
        if len(line) <= 1 or len(line) > 20:
            continue
        else:
            line = line.rstrip()
            print line
        IDs+=1
    result1.set('Total IDs : '+str(IDs))


##

def motion(event):
    x, y = event.x, event.y
    print('{}, {}'.format(x, y))

#root.bind('<Motion>', motion)

              ##



def check():
    letters= ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
    numbers = ["1","2","3","4","5","6","7","8","9","0"]
    minletters = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","k","r","s","t","u","v","x","y","z"]
    data = my_text_box.get("1.0",END).splitlines()

    global theTaxIDs
    global theUniIDs
    global theIprIDs
    global theGOIDs
    global theGeneIDs
    global theKeggIDs
    global thePDBIDs
    global species_list
    global theKeggPaths
    global Words

    theTaxIDs = []
    theUniIDs = []
    theIprIDs = []
    theGOIDs = []
    theGeneIDs = []
    theKeggIDs =[]
    thePDBIDs = []
    species_list = []
    theKeggPaths = []
    theECs = []
    Words = []
    IDs = 0
    for line in data:
        #print line
        line = line.split("\t")[0] # Using first item in line as input
        m_line = line.split(";")
        c_line = line.split(", ")

        if len(line) <= 1:
            pass #
        elif len(line.split(" ")) > 1:
            print "Words : %s" %line
            if line not in species_list:
                species_list.append(line)
                #IDs +=1
                if line not in Words:
                    Words.append(line)
            IDs += 1
        elif line.isdigit() == True: # 5671
            #print " TAX : %s" %line
            if line not in theTaxIDs:
                theTaxIDs.append(line)
                IDs += 1
        elif len(m_line) > 1:
            for line in m_line:
                if line[0] in letters and line[1] in numbers and line[-1] in numbers and line[-2]: # A4HUD2
                    #print " AC : %s" %line
                    if line not in theUniIDs:
                        theUniIDs.append(line)
                        IDs += 1
        elif len(c_line) > 1 : # type(line) is list or type(line) is tuple:
            for lin in line:
                lin = lin[1:-1]
                if lin[0] in letters and lin[1] in numbers and lin[-1] in numbers and lin[-2]: # A4HUD2
                    #print " AC : %s" %lin
                    if lin not in theUniIDs:
                        theUniIDs.append(lin)
                        IDs += 1
        elif line[0] in letters and line[1] in numbers and line[-1] in numbers and line[-2]: # A4HUD2
            # http://www.uniprot.org/help/accession_numbers
            # [OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}
            #print " AC : %s" %line
            if line not in theUniIDs:
                theUniIDs.append(line)
                IDs += 1

        elif line[0:3] == "IPR" and len(line)==9: # IPR002085
            #print " IPR : %s" %line
            if line not in theIprIDs:
                theIprIDs.append(line)
                IDs += 1
        elif line[0:3] == "GO:":# GO:0008270
            #print " GO : %s" %line
            if line not in theGOIDs:
                theGOIDs.append(line)
                IDs += 1
        elif len(line) == 4 and line[0] in numbers and line[1] in letters and line[2] in letters and line[3] in letters : # 1FDH
            #print " PDB : %s" %line
            if line not in thePDBIDs:
                thePDBIDs.append(line)
                IDs += 1
        elif line[0:3] == "EC=":
            #print " Enzyme : %s" %line
            if line not in theECs:
                theECs.append(line)
                IDs +=1
                Words.append(line) # Function FindIN
        #elif line.split("_") >= 1 or line.split(".") >= 1 :
        elif line[0] in letters and line[-1] in numbers and line[-2] in numbers and line[-3] in numbers and line[-4] in numbers: # LINJ_10_0240
            #print " Gene : %s" %line
            if line not in theGeneIDs:
                theGeneIDs.append(line)
                IDs += 1
        elif line[0] in minletters and line[-1] in numbers and len(line) == 8:
            #print " KEGG Path : %s" %line
            if line not in theKeggPaths:
                theKeggPaths.append(line)
                IDs += 1
        #elif len(line)> 0 :
            #print "Word : %s" %line
            

        else:
            print "what is this ID? %s" %line
            if line not in Words:
                Words.append(line)
                IDs += 1

    result2.set('Unique IDs : '+str(IDs))
    print "########################################################################"
    
    print "Unique Tax IDs :" , len(theTaxIDs)
    print "Unique UniProt IDs :" , len(theUniIDs)
    print "Unique InterPro IDs :" , len(theIprIDs)
    print "Unique GO IDs :" , len(theGOIDs)
    print "Unique Gene IDs :" , len(theGeneIDs)
    print "Unique ECs :" , len(theECs)
    print "Unique PDB IDs :" , len(thePDBIDs)
    print "Unique Kegg path IDs :" , len(theKeggPaths)
    print "Unique Words :" , len(Words)
    


def clear_box():
    my_text_box.delete("1.0",END)
    IDs = 0
    result2.set('Unique IDs : '+str(IDs))
    del theInputFiles[:]


# ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt
# UNIPROT HELP
# http://www.uniprot.org/help/accession_numbers



bt11=Button(frame03,text="Check", command=count_input).pack(fill=X) # retrieve_input
#button=Button(frame03,text="Count", command= cipher(text1)).pack()

result1 = StringVar()
result1.set('Total IDS: 0 ')
label1 = Label(frame03,textvariable=result1)
label1.pack()


bt12 =Button(frame03, text=" Ready ", command=check, fg='blue', bd=3, font=('arial',11, 'bold')).pack(fill=X)

result2 = StringVar()
result2.set('Unique IDs : 0 ')
label2=Label(frame03,textvariable=result2)
label2.pack()

bt13 =Button(frame03, text="Clear", command=clear_box).pack(fill=X)


# uni_search_word(theWords)


########   FRAME   #############################33

frame023 = Frame(jp)
frame023.pack()


Label(frame023,text="Search . . . ").pack()

def search_protein_in_uniprot(theUniIDs):
    print "Query Limit = 10 proteins"
    for i in theUniIDs[0:10]:
        webbrowser.open("http://www.uniprot.org/uniprot/"+str(i))




def from_uni_gene_get_kegg_orthologues(theUniIDs,theGeneIDs):

    print "Query Limit = 10 proteins / genes"
    if theUniIDs > 0 :
        for i in theUniIDs[0:10]:
            from_uni_id_get_kegg_orthologues(i)
    elif theGeneIDs > 0:
        for i in theGeneIDs[0:10]:
            from_gene_id_get_kegg_orthologues(i)
    else:
        print "IDs not valid! Try again..."

def from_uni_id_get_kegg_map(theUniIDs):
    print "Query Limit = 10 proteins"
    for i in theUniIDs[:10]:
        from_uni_id_get_kegg_spot(i)
        ##from_uni_ids_get_kegg_path_marker(i)


def from_kegg_path_get_names(theKeggPaths):
    # webbrowser.open("http://www.genome.jp/dbget-bin/get_linkdb?-t+genes+path:"+str(i))
    for i in theKeggPaths:
        webbrowser.open("http://www.genome.jp/dbget-bin/www_bget?pathway+"+str(i))






###  OPEN LINK
bt10 = Button(frame023, text ="Protein UniProt", fg="blue", bd=3 , cursor="hand2", command=lambda: search_protein_in_uniprot(theUniIDs)).pack(side=LEFT) # b1.pack(pady=1)
bt8 = Button(frame023, text="Protein InterPro",  fg="blue", bd=3 , cursor="hand2", command=lambda: from_uni_id_get_interpro_browser(theUniIDs)).pack(side=LEFT) # aux Kegg fucntion
bt9 = Button(frame023, text="Similar Proteins InterPro", fg="blue", bd=3 , cursor="hand2", command=lambda: from_uni_id_find_similar_proteins(theUniIDs)).pack(side=LEFT) # opens WEB
bt7 = Button(frame023, text ="Protein / Gene Kegg orthologues", fg="blue", bd=3 , cursor="hand2", command=lambda: from_uni_gene_get_kegg_orthologues(theUniIDs,theGeneIDs)).pack(side=LEFT) # b1.pack(pady=1)
#bt12 = Button(frame023, text ="Protein Kegg Map", fg="blue", bd=3 , cursor="hand2", command=lambda: from_uni_id_get_kegg_map(theUniIDs)).pack(side=LEFT) # b1.pack(pady=1)
bt12 = Button(frame023, text ="Protein Kegg Map", fg="blue", bd=3 , cursor="hand2", command=lambda: from_uni_ids_get_kegg_path_marker(theUniIDs)).pack(side=LEFT) # b1.pack(pady=1)

bt11 = Button(frame023, text ="Kegg Path Info", fg="blue", bd=3 , cursor="hand2", command=lambda: from_kegg_path_get_names(theKeggPaths)).pack(side=LEFT) # b1.pack(pady=1)
###



# sequence_search_egg

#w = Button(frm1, text="UniProt", command=stats_uniprot_bs, bg="blue", fg="white")
#w.pack(padx=5, pady=5, side=LEFT)



#################################################33
# Same fucntions, just make a for loop

# UNIPROT: theTaxID - > ac_name_ipr_go_path_exi
# INTERPRO: theTaxID -> ac_ipr | ac_ipr_go
# QUICKGO: theUniIDs ->
# (For GuickGO make convert taxid2ipr2go| tax2quickgo| ->
# KEGG: theGeneIDs  -> pathways
# EBI: theIPRIDs -> ipr2go
# BIOBLAST: theUniIDs
# PDB: thePDBID


# UNIQUE TAB FILE:

##################     FRAME  -  RETRIEVE DATA     ################################

frame04 = Frame(jp)
frame04.pack()
goodfont = ('arial', 12)

Label(frame04,text="Retrieve data from Public Databases ...", font = goodfont).pack()

#rooot = Tk()

#height = 5
#width = 5
#for i in range(height): #Rows
    #for j in range(width): #Columns
        #b = Entry(rooot, text="")
        #b.grid(row=i, column=j)
        

##################           CHECKBUTTONS         ############################

databases = ['UniProt', 'BioMart', 'QuickGO', 'KEGG', 'PDB', 'EBI','SwissProt' ]


#  Checkbutton's are associated to several var, kept in list vars

vars = []
for key in databases:
    var = IntVar()
    Checkbutton(frame04, text=key, variable=var, font = goodfont).pack(side=LEFT)
    vars.append(var)

def report():
    """  'Buttons Status'  """
    res = "Database Retrieval :\n"  #princ�pio da msg
    for i in range(len(vars)):
        if vars[i].get() ==1:  # o val�r � 1 (marcado) ou 0 (n�o marcado)
            print databases[i]
            res = res + " " + databases[i] + '\n' #dias � uma lista de strings
    showinfo("Escolhas", res)




def convert_ac2gene(theTaxID):
    #access_uniprot_with_biomart(theTaxID)

    theFilename = "5671_ac_entry_ipr_name_gene_go_path_local_exi.txt"
    """Organism ID	Entry name	Entry	InterPro	Protein names	Gene names	Gene ontology IDs	Gene ontology (GO)	Pathway	Protein existence	Organism
    5671
    ASNA_LEIIN
    A4HUY0
    IPR025723; IPR016300; IPR027542; IPR027417;
    ATPase ASNA1 homolog (EC 3.6.-.-) (Arsenical pump-driving ATPase homolog) (Arsenite-stimulated ATPase)
    LinJ11.0710 LinJ_11_0710
    GO:0005524; GO:0016887; GO:0005783; GO:0046872; GO:0045048; GO:0006810
    ATPase activity; ATP binding; endoplasmic reticulum; metal ion binding; protein insertion into ER membrane; transport		Inferred from homology	Leishmania infantum"""
    f = open(theFilename, "r")
    for i in f.readlines():
        #print i
        acID = i.strip().split("\t")[2]
        geneIDs = i.strip().split("\t")[5]
        print acID
        for j in geneIDs.split(" "):
            if "Linj" in j:
                print j
            elif "LINJ" in j:
                print j
            elif "LinJ" in j:
                print j
            elif "" == j:
                print "NO gene NAME", acID
            else:
                print "Gene name:", j

def lab2uni2ipr2go2kegg(theFilename):
    f=open(theFilename, "r")
    for i in f.readlines():
        print i




"""S = Scrollbar(jp)
    T = Text(jp, height=4, width=50)
    S.pack(side=RIGHT, fill=Y)
    T.pack(side=LEFT, fill=Y)
    S.config(command=T.yview)
    T.config(yscrollcommand=S.set)
    #T.insert(END, quote)"""
# TypeError: search() got an unexpected keyword argument 'format'

from run_bioblast_tk import *

def print_data():



    """Start data retrieval by check button """
    #theTaxID = E01.get()
    #print "Tax ID is ", theTaxID
    print "###################################      RUNNING      ###################################"
    print theTaxIDs
    print theUniIDs
    print theIprIDs
    print theGOIDs
    print theGeneIDs
    print theKeggIDs
    print thePDBIDs

    #theGeneIDs = convert_ac2gene(the_aclist)

    #theUniProtIDs =
    # define 3 last functions glue
    for i in range(len(vars)):
        if vars[i].get() == 1:
            print "\n > > > Database Loading . . .", databases[i]
            if databases[i] =="UniProt":
                if len(theTaxIDs) > 0 :
                    for id in theTaxIDs:
                        from_tax_id_get_bs_uniprot(id) 
                
                elif len(theUniIDs) < 500 and len(theUniIDs) >0: 
                    expected = 0.0036
                    minutes = len(theUniIDs)*float(expected)
                    m = float("{0:.1f}".format(minutes))
                    hours = minutes /60
                    h = float("{0:.1f}".format(hours))
                    response = "Time Expected: " , m, " minutes,  ",h, " hours."
                    print response
                    from_uni_id_get_bs_uniprot1(theUniIDs)

                elif len(theUniIDs) > 499:
                    expected = 0.0036
                    minutes = len(theUniIDs)*float(expected)
                    m = float("{0:.1f}".format(minutes))
                    hours = minutes /60
                    h = float("{0:.1f}".format(hours))
                    response = "Time Expected: " , m, " minutes,  ",h, " hours."
                    print response

                    Freq = 1000
                    Dur = 100
                    winsound.Beep(Freq,Dur)
                    quest = raw_input("Proceed ? Y / N : ")

                    choicesy = ["y", "Y"]
                    choicesn = ["N", "n"]
                    if quest in choicesn:
                        print "Try another database ! Thank You "

                    elif quest in choicesy:
                        from_uni_id_get_bs_uniprot1(theUniIDs)

                #elif len(theUniIDs) > 0:
                    #from_uni_id_get_bs_uniprot1(theUniIDs) # RETRIVE ONE FILE
                elif len(theGeneIDs) > 0:
                    from_gene_id_get_bs_uniprot1(theGeneIDs)
                    
                else:
                    print "IDs are not valid  for UniProt. . ."

            if databases[i] =="BioMart":
                if len(theTaxIDs) > 0 :
                    for id in theTaxIDs:
                        from_tax_id_get_bs_interpro2go(id)
                elif len(theUniIDs) < 300 : # 0 and < 250:
                    from_uni_id_get_bs_ac2ipr2go(theUniIDs)
                elif len(theUniIDs) > 299:
                    expected = 0.052
                    minutes = len(theUniIDs)*float(expected)
                    m = float("{0:.1f}".format(minutes))
                    hours = minutes /60
                    h = float("{0:.1f}".format(hours))
                    response = "Time Expected: " , m, " minutes,  ",h, " hours."
                    print response
                    Freq = 1000
                    Dur = 100
                    winsound.Beep(Freq,Dur)
                    quest = raw_input("Proceed ? Y / N : ")
                    choicesy = ["y", "Y"]
                    choicesn = ["N", "n"]
                    if quest in choicesn:
                        print "Try another database ! Thank You "

                    elif quest in choicesy:
                        from_uni_id_get_bs_ac2ipr2go(theUniIDs)


                elif len(theIprIDs) > 0 :
                    from_ipr_id_get_bs_ipr2go(theIprIDs)
                else:
                    print "IDs are not valid for BioMart. . ."

            if databases[i] =='QuickGO':
                if len(theTaxIDs) > 0:
                    for id in theTaxIDs:
                        from_tax_id_get_quickgo(id)
                elif len(theUniIDs) > 0:
                    from_uni_ids_get_quickgo(theUniIDs)
                elif len(theGOIDs) > 0 :
                    from_go_ids_get_quickgo(theGOIDs)
                else:
                    print "IDs are not valid  for QuickGO. . ."

            if databases[i] =='KEGG':

                if len(theUniIDs) > 0:
                    from_uni_ids_get_kegg_paths(theUniIDs)
                    #from_uni_ids_get_kegg(theUniIDs)
                    #for i in theUniIDs:
                    #from_uni_ids_get_kegg_path(i)
                #elif len(theUniIDs) < 6:
                    #for i in theUniIDs:
                        #from_uni_id_get_kegg(i)

                elif len(theGeneIDs) > 0: # theKeggIDs
                    if len (theGeneIDs) < 20:
                        #from_gene_uni_ids_get_kegg_image(i)
                        #else:
                        for id in theGeneIDs:
                            from_gene_ids_get_kegg_path(id)
                            ##from_kegg_ids_get_kegg_path(id)
                            #theKeggIDs.join("\n")
                            #access_kegg_with_biomart(theGeneIDs)#bioservices_functions_tk-
                            #convert_ac2gene(theTaxID)

                else:
                    print "IDs are not valid for KEGG. . ."

            if databases[i] =='PDB':
                if len(thePDBIDs) > 0 :
                    for id in thePDBIDs:
                        from_pdb_id_get_pdb(id)
                elif len(theUniIDs) > 0 :
                    for id in theUniIDs:
                        getPdb(id)
                    #from_uni_ids_get_pdb(theUniIDs)
                    getPdbinfo(theUniIDs)
                else:
                    print "IDs are not valid for PDB. . ."

            if databases[i] =='SwissProt':
                if len(theUniIDs) > 0 :
                    expected = 0.008
                    minutes = len(theUniIDs)*float(expected)
                    m = float("{0:.1f}".format(minutes))
                    hours = minutes /60
                    h = float("{0:.1f}".format(hours))
                    response = "Time Expected: " , m, " minutes,  ",h, " hours."
                    print response
                    get_bioblast_record(theUniIDs)
                    # access_bioblast_with_biopython(theUniIDs) # TRANFER TO SUMMARY Analysis as BioBlas

            else:
                print "No Databases selected ." 
    print "######################################     Complete!     ################################################"


                        # import run_biopython_tk
                        #access_blast_with_biopython(theTaxID)


import run_biopython_tk
# uni_search_words(theTerms)

# The  most important value to compare between databases is also the most present and makes it easy to find.
# (HELP UNIPROT LIST)

#Button(jp, text='Estado dos Bot�es', command=report, font = goodfont).pack(fill=X) # save data
# Button(jp, text='Sair', command=exit).pack(fill=X) ERROR




def file_save_x():
    f = tkFileDialog.asksaveasfile(mode='w', defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        return
    text2save = str(text.get(1.0, END)) # starts from `1.0`, not `0.0`
    for i in text2save:
        f.write(i)
    f.close() # `()` was missing.
    filename3= "" #taxid2ipr2go - action merged
    

def save_bmipr_taxid2ipr2go():
    f = tkFileDialog.asksaveasfile(mode='w', defaultextension=".txt")
    if f is None: # asksaveasfile return `None` if dialog closed with "cancel".
        showinfo("Warning", "No File Saved.")
    #text2save = str(text.get(1.0, END)) # starts from `1.0`, not `0.0`
    #f.write(text2save)
    #f.close() # `()` was missing.
    filename3= f # taxid2ipr2go - action merged
    access_bmipr_ac2IPR2go(filename1,filename2,filename3)


###########################################################################


bt7 = Button(text ="<<<    Go    >>>", command=print_data, fg="blue", bd=3, font=('arial',11, 'bold')).pack(fill=X) # b1.pack(pady=1)

#theInputFiles = ["output/2pep_260_filter_260_tax5671_bs_uniprot.txt","output/2pep_246_down.txt", "output/2pep_246_down.txt"]

###################################################################################3
#theInputFiles = ["5671_ac_entry_ipr_name_gene_go_path_local_exi.txt","5671_bmq_ac_IPR_go.txt","qui_uni_goid_type_db.txt"]
def save_unique(thefiles):
    files = ["theInputList1", "theInputList2", "theInputList3", "theInputList4","theInputList5","theInputList6","theInputList7","theOutput"]
    vars = []
    for key in files:
        #var = IntVar()
        vars.append(key)
    count = 1
    #print theFilenames1
    for i in range(len(thefiles[:3])):
        if len(thefiles[i]) > 1:
            vars[i] = thefiles[i]
            print i, vars[i]
        count += 1
    theOutput = str(raw_input("Save file as : ")+"_merged.txt")
    save_input_file_to_append_by_ac_ids(vars[0], vars[1], vars[2],vars[3],vars[4],vars[5],vars[6], theOutput)
    #del theInputFiles[:]

#save_unique(theInputFiles)

###############33





# bt5 = Button(text =" Merge by AC ID ", command=save_unique,font=('arial',11, 'bold')).pack(fill=X) # check ok button

#bt5 = Button(text ="Save Unique tab File ... ", command=file_save).pack(fill=X) # check ok button



"""global _run_once
global refresh
_run_once=0

def choi():
    print "helo"
def refresh_clicked(event):
    global _run_once
    global refresh
    #def wrapper():
    #if not wrapper.has_run:
    #wrapper.has_run = True
    #"_run_once =0
    if _run_once == 0:
        #framex = Frame(framex, borderwidth=5)
        #framex.grid()
        bt20= tk.Button( text="LAcFilter", command= choi , fg="black") #bt20.pack( padx=5, pady=1, side=LEFT)
        bt20.grid(row=0, column = 0,padx=1, pady=1)

        #refresh_clicked.func_code = (lambda:None) 
        _run_once = _run_once + 1
        refresh.destroy()

        



#framex = Frame(jp)
#framex.pack()



#et1 = tk.Entry(frame)
#et1.insert(0, 10)
#et1.grid(row=0,column=0,sticky=tk.W)
#label_contents = tk.StringVar()
#label_contents.set(et1.get())

        
        
#tk.Label(frame, textvariable=label_contents).grid(row=1, column=0, sticky=tk.W)
#refresh = tk.Button(framex, text='Refresh', command = refresh_clicked)
#refresh.grid(row=2, column=0, sticky=tk.W)
#refresh.destroy()

#per comments:
"""


##############


##################################################3
##### summary_analysis
############################################3#####

#def raise_frame(frame):
    #frame.grid()
    #bt= Button(text ="Summ ",command=faznada,bd=2) # check ok button
    #bt.pack()
    
    #frame.tkraise()


#if1 = Frame(jp)


#for frami in (if1, if2, if3, if4):
    
    #frami.grid(row=0, column=0, sticky='news')


#lb5 = Button(text ="Summary Analysis ",command=lambda: raise_frame(if1),bd=2) # check ok button
lb5 = Label(text ="Summary Analysis ",bd=2) # check ok button

lb5.pack()


###COLOURS
"""For example, "#fff" is white, "#000000" is black, "#000fff000" is pure green, and "#00ffff" is pure cyan (green plus blue).

You can also use any locally defined standard color name. The colors "white", "black", "red", "green", "blue", "cyan", "yellow", and "magenta" will always be available."""


####
theFilenames = ["outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt",
                "outputs/5671_bmq_ac_ipr_go.txt",
                "",
                "",
                "",
                ""]

def summary_funtions_1(theFilenames):

    print "\nUniProt by tax id "
    print '[1] Loading protein accessions '
    print theFilenames[0]
    Lin_Proteins = load_proteins_accessions(str(theFilenames[0])) # "output/5671_ac_entry_ipr_name_gene_go_path_local_exi_short.txt")

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
    print '[3] Map InterPro 2 GO with BiomartView '


####





def summary_analysis():
    "Loading data . . . . "
    i = 0               # usually has a control variable to start the while loop
    while i < 20:       # terminating condition: while loop will terminate if i >= 10
        i += 1          # changing the value of the control variable
        print '|' * i
    print "                                                        100%"

    import run_summary_functions_bs_tk  # loads and runs


# run_uniprot_summary_functions_bs()
# run_interpro_summary_functions_bs()
# run_quickgo_summary_functions_bs()
# run_kegg_summary_functions_bs()
#


#from run_summary_functions_manual_tk import *
# FIND: load_ipr2go_list'

def run_res():
    run_interpro_summary_functions_bs()
    run_uniprot_summary_functions_bs()




#from raio_rotacao_tk import* # raio?

#b1 = Button(jp, text ="Summary Analysis", command=run_uniprot_summary_functions_bs()).pack(fill=X) # check ok button

#b3 = Button(jp, text ="Summary Analysis 2 ", command=run_res).pack(fill=X) # check ok button


###############

from run_summary_functions_manual_tk import * #stats_bmv_manual
#from run_summary_functions_bs_tk import * # stats_uniprot_bs, stats_interpro_bs

import summary_functions_tk

frm1 = Frame(jp)
frm1.pack()

# import summary_functions_tk
def retrieve2stats(theInputFiles):
    for i in theInputFiles:
        print 'Reading tab file -----------------------------------------------'
        filename = i
        prot_data, fields = readUniProtTab_sum(filename, with_fields=True)
        print 'there are', len(prot_data), 'proteins in file', filename



def stats_uniprot_bs():
    summary_funtions_1(theFilenames)


#w = Button(frm1, text=" UniProt | InterPro | QuickGO | Kegg | BioBlast ", command= lambda: retrieve2stats(theInputFiles), bd=5, bg="blue", fg="white") # command=stats_uniprot_bs
#w.pack(padx=5, pady=5, side=LEFT)

# TheInputFiles if used in the MERGED FUNCTION will not work for RETRIVE2STATS
# For that reason SUMMARY will be in another frame, selecting one file at the time!

################################################################################################################
########################################## SUM FRAME  ##########################################################
################################################################################

frm22 = Frame(jp)
frm22.pack()



theFilenames1 = [] #["_bmq_ipr_go.txt", "outputs/5671_ac_entry_ipr_name_gene_go_path_local_exi.txt", "outputs/5671_bmq_ac_ipr_go.txt"]


##############################################################


# MultiListBox - View Mode


# http://code.activestate.com/recipes/52266/
class MultiListbox(Frame):
    def __init__(self, master, lists):
        Frame.__init__(self, master)
        self.lists = []
        for l,w in lists:
            frame = Frame(self); frame.pack(side=LEFT, expand=YES, fill=BOTH)
            Label(frame, text=l, borderwidth=1, relief=RAISED).pack(fill=X)
            lb = Listbox(frame, width=w, borderwidth=0, selectborderwidth=0,

            relief=FLAT, exportselection=FALSE)
            lb.pack(expand=YES, fill=BOTH)
            self.lists.append(lb)
            lb.bind('<B1-Motion>', lambda e, s=self: s._select(e.y))
            lb.bind('<Button-1>', lambda e, s=self: s._select(e.y))
            lb.bind('<Leave>', lambda e: 'break')
            lb.bind('<B2-Motion>', lambda e, s=self: s._b2motion(e.x, e.y))
            lb.bind('<Button-2>', lambda e, s=self: s._button2(e.x, e.y))

        frame = Frame(self); frame.pack(side=LEFT, fill=Y)
        Label(frame, borderwidth=1, relief=RAISED).pack(fill=X)
        sb = Scrollbar(frame, orient=VERTICAL, command=self._scroll)
        sb.pack(expand=YES, fill=Y)
        bb = Scrollbar(frame,orient=HORIZONTAL, command = self._scroll)
        bb.pack(side=BOTTOM, expand=YES, fill=X)
        self.lists[0]['yscrollcommand']= sb.set
        self.lists[0]['xscrollcommand'] = bb.set

    def _select(self, y):
        row = self.lists[0].nearest(y)
        self.selection_clear(0, END)
        self.selection_set(row)
        return 'break'

    def _button2(self, x, y):
        for l in self.lists: l.scan_mark(x, y)
        return 'break'

    def _b2motion(self, x, y):
        for l in self.lists: l.scan_dragto(x, y)
        return 'break'

    def _scroll(self, *args):
        for l in self.lists:
            apply(l.yview, args)
        return 'break'

    def curselection(self):
        return self.lists[0].curselection()

    def delete(self, first, last=None):
        for l in self.lists:
            l.delete(first, last)

    def get(self, first, last=None):
        result = []
        for l in self.lists:
            result.append(l.get(first,last))

        if last: return apply(map, [None] + result)
        return result

    def index(self, index):
        self.lists[0].index(index)

    def insert(self, index, *elements):
        for e in elements:
            i = 0
            for l in self.lists:
                l.insert(index, e[i])
                i = i + 1
                

    def size(self):
        return self.lists[0].size()

    def see(self, index):
        for l in self.lists:
            l.see(index)

    def selection_anchor(self, index):
        for l in self.lists:
            l.selection_anchor(index)

    def selection_clear(self, first, last=None):
        for l in self.lists:
            l.selection_clear(first, last)

    def selection_includes(self, index):
        return self.lists[0].selection_includes(index)

    def selection_set(self, first, last=None):
        for l in self.lists:
            l.selection_set(first, last)



def preview(aFile):

    tk = Tk()
    tk.title = "View Mode"
    toplabel = "Filename : "+ str(aFile)
    Label(tk, text=toplabel).pack()
    theFile = open(aFile, "r")
    header = theFile.readline()
    header = header.split("\t") #print header
    mlb = MultiListbox(tk, ((str(header[0]), 10), (str(header[1]), 10), (str(header[2]), 10),(str(header[3]), 40),(str(header[4]), 10),
                            (str(header[5]), 10),(str(header[6]), 10),(str(header[7]), 5),(str(header[8]), 5),(str(header[9]), 5)))#,
                            #(str(header[15]), 10),(str(header[16]), 10),(str(header[17]), 10),(str(header[18]), 10)))
    #mlb = MultiListbox(tk, ((str(header[0]), 10), (str(header[1]), 10), (str(header[2]), 10),(str(header[3]), 10),(str(header[4]), 10),(str(header[5]), 10),(str(header[6]), 10)))
    #for i in header:
        #mlb = MultiListbox(tk, (str(i),10))
    #mlb = MultiListbox(tk, (('Subject', 10), ('Sender', 10), ('Date', 10)))
    #for i in range(1000):
    for i in theFile.readlines()[1:]:
        line = i.strip()
        #print line 
        items = line.split("\t")
        
        if len (items) >2 :#and len(items)<15:
            #print line
            #print items
            #for j in len(items):#mlb.insert(END, ('Important Message: %d' % i, 'John Doe', '10/10/%04d' % (1900+i)))
            mlb.insert(END, items)

    mlb.pack(expand=YES,fill=BOTH)
    tk.mainloop()


#preview("all_kk.txt") 


##########################################
#de#f poll():
  #  lab.after(200, poll)
  #  sel = L.curselection()
  #  lab.config(text=str(sel))

#import tkinter

def close_window(new_window):
    new_window.destroy()

def open_col_list():

    import Tkinter
    F1 = Tkinter.Frame()
    s = Tkinter.Scrollbar(F1)
    L = Tkinter.Listbox(F1)

    s.pack(side=Tkinter.RIGHT, fill=Tkinter.Y)
    L.pack(side=Tkinter.LEFT, fill=Tkinter.Y)

    s['command'] = L.yview
    L['yscrollcommand'] = s.set

    for i in range(30):
        L.insert(Tkinter.END, str(i))



    #poll()


    F1.pack(side=Tkinter.TOP)

    F2 = Tkinter.Frame()
    lab = Tkinter.Label(F2)

    lab.pack()
    F2.pack(side=Tkinter.TOP)

    lab.after(200, poll)
    sel = L.curselection()
    lab.config(text=str(sel))

    Tkinter.mainloop()



###################################################################################################

# http://effbot.org/tkinterbook/scrollbar.htm

def CurSelet(evt):

    #values = [mylistbox.get(idx) for idx in mylistbox.curselection()]
    #print ', '.join(values)
    numbers = ["1", "2","3","4","5","6","7","8","9","0"]
    value=str(mylistbox.get(mylistbox.curselection()))#mylistbox.get(ACTIVE)))

    print value
    items = value.split("    ")
    #print items[1]
    #for i in items:
    #   print i
    if items[1][:3] == "IPR":
        webbrowser.open("http://www.ebi.ac.uk/interpro/entry/"+str(items[1]))
    elif items[1][1] in numbers: # and items[1][0] in letters
        webbrowser.open("http://www.ebi.ac.uk/interpro/protein/"+str(items[1]))
    elif items[1][:3] == "GO:":
        webbrowser.open("http://www.ebi.ac.uk/interpro/search?q="+str(items[1]))


#################################################  LIST BOX  ##########################################

def open_list_box1(filename):

    #from tkinter import *
    #from tkinter import ttk
    # http://stackoverflow.com/questions/7616541/get-selected-item-in-listbox-and-call-another-function-storing-the-selected-for

    root1=Tk()
    #new_window = Toplevel(root)
    root1.title(str(filename))
    sizex = 1200
    sizey = 400
    posx  = 90
    posy  = 80
    root1.wm_geometry("%dx%d+%d+%d" % (sizex, sizey, posx, posy))
    #itemsforlistbox=['one','two','three','four','five','six','seven']
    #items = "    ".join(itemsforlistbox)

    scrollbar = Scrollbar(root1)
    scrollbar.pack(side=RIGHT, fill=Y)

    global mylistbox
    mylistbox=Listbox(root1,width=250,height=80,font=('times',13), yscrollcommand=scrollbar.set)
    mylistbox.place(x=32,y=90)
    mylistbox.bind("<Double-Button-1>", CurSelet) # '<<ListboxSelect>>'


    #button = Button(root1, text='push') #, command=get_item(value))
    #button.pack(side="CENTER")


    #for it in items:
    choic = open(filename,"r")
    count = 0
    for i in choic.readlines():
        #print i
        i = i.strip()
        items = i.split("\t")
        new = "    ".join(items)
        neww = str(count) + "    " + new
        mylistbox.insert(END,neww)
        count += 1



    mylistbox.pack(side=LEFT, fill=BOTH)


    scrollbar.config(command=mylistbox.yview)
    choic.close()

    #ok= Button(root1, text="Close",command=lambda: close_window(root1)).pack(side=BOTTOM)


    root1.mainloop()



################################################     COMBOBOX    ###################################################3

def clean_combo():
    del theFilenames1[:]
    del theInputFiles[:]
    box_value = StringVar()
    cbox.set("")
    box_value2 = StringVar()
    cbox2.set("")
    #cbox['values'] = files
    #cbox.bind("<<ComboboxSelected>>", newselection)

def clean_combo2():
    del theFilenames1[:]
    del theInputFiles[:]
    box_value2 = StringVar()
    cbox2.set("")




# http://stackoverflow.com/questions/18906047/how-to-update-values-to-the-listbox-under-combobox-in-ttk-python33
#update list upon drop down
files = ()
global box_value
global box_value2
#choice = ''



def updatecombolist(theFilenames1):

    for i in theInputFiles:
        if i not in theFilenames1:
            theFilenames1.append(i)

    files = tuple(theFilenames1)
    global choice
    #global box_value
    #global theFilenames1
    #list = self.getPortLst()

    box_value = StringVar()

    #cbox = ttk.Combobox(frm22, textvariable=box_value)
    if len(theFilenames1) == 0:
        choice = "null"
    else:
        choice = theFilenames1[-1]
        #print cbox.set(choice)
        #print cbox.get()
        #choice = cbox.get()#print choice = choice
        #print choice
        #cbox
        #print box_value.get()
        box_value = StringVar()
        box_value = choice
        #print "Seleced : ", choice,"\n"
        cbox.set(choice)
        #print box_value # PY_VAR19, PY_VAR28

    cbox['values'] = files
    cbox.bind("<<ComboboxSelected>>", newselection)



    #box.value.set()
    #textvariable.set(box_value)
    #print textvariable # _name_ac1.txt


    cbox.current()


    #choice = cbox.get()
    #print choice
    #print box_value # PY_VAR16
    #print type(choice)
    #box.grid(column=0, row=0)
    #cbox.pack(side=TOP)



def updatecombolist2(theFilenames1):

    for i in theInputFiles:
        if i not in theFilenames1:
            theFilenames1.append(i)

    files = tuple(theFilenames1)
    global choice2
    #global box_value
    #global theFilenames1
    #list = self.getPortLst()

    box_value2 = StringVar()

    #cbox = ttk.Combobox(frm22, textvariable=box_value)
    if len(theFilenames1) == 0:
        choice2 = "null"
    else:
        choice2 = theFilenames1[-1]
        #print cbox.set(choice)
        #print cbox.get()
        #choice = cbox.get()#print choice = choice
        #print choice
        #cbox
        #print box_value.get()
        box_value2 = StringVar()
        box_value2 = choice2
        #print "Seleced : ", choice2,"\n"
        cbox2.set(choice2)
        #print box_value # PY_VAR19, PY_VAR28

    cbox2['values'] = files
    cbox2.bind("<<ComboboxSelected>>", newselection2)



    #box.value.set()
    #textvariable.set(box_value)
    #print textvariable # _name_ac1.txt


    cbox2.current()


    #choice = cbox.get()
    #print choice
    #print box_value # PY_VAR16
    #print type(choice)
    #box.grid(column=0, row=0)
    #cbox.pack(side=TOP)




##################################################################################################
###########################################     COMBOBOX 0    #####################################

thefile = ""
def delete1(thefile):

    box_value = StringVar()
    cbox.set("")

    if thefile in theFilenames1:
        theFilenames1.remove(thefile)
        if thefile in theInputFiles:
            theInputFiles.remove(thefile)
    print thefile, " is removed !"

#

def delete2(thefile):


    box_value2 = StringVar()
    cbox2.set("")
    if thefile in theFilenames1:
        theFilenames1.remove(thefile)
        if thefile in theInputFiles:
            theInputFiles.remove(thefile)
    print thefile, " is removed !"
    #else:
    #print thefile, " is removed !"

#updatecombolist2(theFilenames1)

#####################################################################################################################

def paste_text(choice):

    aFile = open(choice)
    my_text_box.delete(1.0 , END)
    for line in aFile.readlines():
        #print line
        my_text_box.insert( END, line +"\n")
    aFile.close()

######################################################################################################################
###################################################   FRAME    ##############################################################

frm222 = Frame(jp)
frm222.pack()

label33 = ttk.Label(frm222, text=' Select Background :')
label33.pack(side=LEFT)


box_value2 = StringVar()
#box_value.set("No File Selected!")
cbox2 = ttk.Combobox(frm222, width = 65, textvariable=box_value2, postcommand = lambda: updatecombolist2(theFilenames1))#,state='readonly')
cbox2.pack(side=LEFT,anchor=N,expand=YES)


w2 = Button(frm222, text=" Upload File ", command= upload_file_list2, bd=2 )# , bg="blue", fg="white")
w2.pack(padx=5, pady=5, side=LEFT)

b000= Button(frm222, text="Preview", command= lambda: preview(choice2), bd=2).pack(side=LEFT)# lambda: open_text_list(choice2),
b010= Button(frm222, text="Edit", command= lambda: open_file_edit(choice2), bd=2).pack(side=LEFT)# lambda: open_text_list(choice2)
#opennotepadFile


b004= Button(frm222, text="Del", command= lambda: delete2(choice2), bd=2).pack(side=LEFT) # lambda: open_text_list(choice)



# http://www.tutorialspoint.com/python/tk_listbox.htm
#b01 = Button(frm22, text ="TAB", command = open_col_list, bd=2).pack(side=LEFT)#open_list_box

w03 = Button(frm222, text="Clear", command= clean_combo)
w03.pack(side=RIGHT)



#####
###########################################     COMBOBOX 1      #####################################
# http://www.tkdocs.com/tutorial/widgets.html
# https://mail.python.org/pipermail/tkinter-discuss/2010-April/002214.html


label222 = ttk.Label(frm22, text=' Select File :')
label222.pack(side=LEFT)
box_value = StringVar()
#box_value.set("No File Selected!")
cbox = ttk.Combobox(frm22, width = 65, textvariable=box_value, postcommand = lambda: updatecombolist(theFilenames1))#,state='readonly')
cbox.pack(side=LEFT)


w2 = Button(frm22, text=" Upload File ", command= upload_file_list, bd=2 )# , bg="blue", fg="white")
w2.pack(padx=5, pady=5, side=LEFT)

b00= Button(frm22, text="View", command= lambda: open_list_box1(choice), bd=2).pack(side=LEFT) # lambda: open_text_list(choice),open_list_box1(choice)

b005= Button(frm22, text="Paste", command= lambda: paste_text(choice), bd=2).pack(side=LEFT) # lambda: open_text_list(choice)
b0155= Button(frm22, text="Edit", command= lambda: open_file_edit(choice), bd=2).pack(side=LEFT)# lambda: open_text_list(choice2)

###btt8 = Button(frm22, text="Open", command= lambda: open_excel_macro(choice), bd=2).pack(side=LEFT) # lambda: open_text_list(choice)



b003= Button(frm22, text="Del", command= lambda: delete1(choice), bd=2).pack(side=LEFT) # lambda: open_text_list(choice)

# http://www.tutorialspoint.com/python/tk_listbox.htm
#b01 = Button(frm22, text ="TAB", command = open_col_list, bd=2).pack(side=LEFT)#open_list_box

w3 = Button(frm22, text="Clear", command= clean_combo)
w3.pack(side=RIGHT)

#######################3########################################################################################

# ************** LISTBOX ****************

def onselect(evt):
    # Note here that Tkinter passes an event object to onselect()
    w = evt.widget
    index = int(w.curselection()[0])
    value = w.get(index)
    print 'You selected item %d: "%s"' % (index, value)

def open_list_box():

    root3 = Tk()
    lb = Listbox(root3, name='lb')
    lb.bind('<<ListboxSelect>>', onselect)





############################################################################################
    
LINKS=("http://www.python.org", "http://www.heaven.com")

def showLink(event):
    idx= int(event.widget.tag_names(CURRENT)[1])
    print LINKS[idx]

def open_text_list(choice):

    root2 = Tk()
    # https://mail.python.org/pipermail/tkinter-discuss/2008-August/001618.html
    global selected
    txt=Listbox(root2)

    #txt.insert(END, "Press ")
    #txt.insert(END, "here ", ('link', str(0)))
    #txt.insert(END, "for Python. Press ")
    #txt.insert(END, "here ", ('link', str(1)))
    #txt.insert(END, "for Heaven.")
    #txt.tag_config('link', foreground="blue")
    #txt.tag_bind('link', '<Button-1>', showLink) # use selectd

    with open(choice,"r") as f:
        for i in f.readlines():
            #print i
            #i = i.split("\t")
            txt.insert(END, i)
    selected = txt.curselection()
    #print selected
    f.close()
    txt.pack(expand=True, fill="both")
    root2.mainloop()


def use_selected(selected):

    print selected
    for i in selected.split("\t"):
        webbrowser.open("http://www.ebi.ac.uk/interpro/protein/"+str(i)) # "http://www.python.org") #+str(selected))
        #from_uni_id_get_interpro_browser(theUniIDs):
        #print "Query Limit = 10 proteins "
        #for i in theUniIDs[0:10]:
        #   webbrowser.open("http://www.ebi.ac.uk/interpro/protein/"+str(i))





def clear_sum():
    del theFilenames1[:]



def clickAbout():
    toplevel = Toplevel()
    label1 = Label(toplevel, text=ABOUT_TEXT, height=0, width=100)
    label1.pack()
    label2 = Label(toplevel, text=DISCLAIMER, height=0, width=100)
    label2.pack()

################################################################################################3

# helps : http://www.tkdocs.com/widgets/combobox.html

def combo1(frm2,theFilenames1):
    value = StringVar()
    box = ttk.Combobox(frm2, textvariable=value, state='readonly')
    files = tuple(theFilenames1)
    box['values'] = files
    box.current(0)
    #box.grid(column=0, row=0)
    box.pack(side=LEFT,expand=YES)

##############################################################################################3

def combo2(frm222,theFilenames1):
    value = StringVar()
    box = ttk.Combobox(frm2, textvariable=value, state='readonly')
    files = tuple(theFilenames1)
    box['values'] = files
    box.current(0)
    #box.grid(column=0, row=0)
    box.pack(side=LEFT,expand=YES)

###
# http://stackoverflow.com/questions/6876518/set-a-default-value-for-a-ttk-combobox
# http://stackoverflow.com/questions/17757451/simple-ttk-combobox-demo
#class App:
#def __init__(self, parent):
 #   self.parent = parent
  #  self.value_of_combo = 'X'
   # self.combo()

def newselection( event):
    value_of_combo = cbox.get()
    print "#################################################################"
    print(value_of_combo), " #  >>>>>    New Selection (1) !"
    global choice
    choice = cbox.get()
    #print choice
    #cbox = ttk.Combobox(frm2, textvariable=choice, state='readonly')
    #box_value = choice#cbox.get()
    #cbox.pack(side=TOP)

def newselection2( event):
    value_of_combo = cbox2.get()
    print "#################################################################"
    print(value_of_combo), " #  >>>>>    New Selection (2) !"
    global choice2
    choice2 = cbox2.get()
    #print choice
    #cbox = ttk.Combobox(frm2, textvariable=choice, state='readonly')
    #box_value = choice#cbox.get()
    #cbox.pack(side=TOP)


def comboxi():
    box_value = StringVar()
    box = ttk.Combobox(frm22, textvariable=box_value)
    box.bind("<<ComboboxSelected>>", newselection)
    # ...
###
def from_file_get_stats(choice): # Counter @ summary functions

    #for i in theFilenames1:
    #print str(choice) # _name_ac1.txt
    #print type(choice) # unicode
    print '\nReading tab file -----------------------------------------------'
    filename = choice
    prot_data, fields = readUniProtTab_sum(filename,n=0, with_fields=True)
    print '\nThere are', len(prot_data), 'proteins in file', filename,"\n#####################################################\n"



#frm2 = Frame(jp)
#frm2.pack()

#status = Label(jp, text="Preparing to do nothing...", bd=1, relief=SUNKEN, anchor=W)
#status.pack(side=BOTTOM,fill=X)


#photo = PhotoImage(file="icon.gif")
#w = Label(parent, image=photo)
#w.photo = photo
#w.pack()

#w3 = Label(frm2, text = print_filenames)
#w3.pack(side=LEFT)

#photo = PhotoImage(file="icon.gif")
#w = Label(parent, image=photo)
#w.photo = photo
#w.pack()

"PAth interactions by IPR" # http://www.ebi.ac.uk/interpro/entry/IPR000796/pathways-interactions



############################################################################################################

###            ********************************* SUMMARY   *************************************************
###                    **********************    BUTTONS   *****************************


goodfont = ('arial', 12)

frm8 = Frame(jp)
frm8.pack()#frm.pack(side=LEFT)



#Label(frm1,text="Get Unique Combined Data File in ...", font = goodfont).pack(fill=X)

#theFilenames1

def get_gene_lists():

    for i in theInputFiles:
        print i
def get_interpro2go_domains():

    # http://sebastianraschka.com/Articles/2014_scikit_dataprocessing.html

    import csv
    import urllib

    url = 'https://raw.githubusercontent.com/rasbt/pattern_classification/master/data/wine_data.csv'
    csv_cont = urllib.request.urlopen(url)
    csv_cont = csv_cont.read() #.decode('utf-8')

    # Optional: saving the data to your local drive
    with open('./wine_data.csv', 'wb') as out:
        out.write(csv_cont)



# ***************** IMAGE **********************

#photo = PhotoImage(file="ewds.png")
#label = Lable(jp, image=photo)
#label.pack()

#b0 = Button(frm1, text ="Table", command= lambda: open_text_list(choice), bg="blue", fg="white").pack(padx=5, pady=5, side=LEFT)# ()print_data

def upload_bg():

    filename_input = tkFileDialog.askopenfilename(title="Open File Background", filetypes=[("File TXT",".txt"),("All files",".*")])
    if len(filename_input) > 0:
        showinfo("OK", "Filename :\n"+filename_input)
    else:
        showinfo("Warning", "No file was selected.")
    #aFile = open(filename_input)
    #my_text_box.delete(1.0 , END)
    #for line in aFile.readlines():
        #print line
        #my_text_box.insert( END, line +"\n")
    name = os.path.basename(filename_input)
    bg_choice = name

    return bg_choice

    #enrichment_analysis_domains

####################################################################
####################################################################

#newf = from_tax_id_get_bs_interpro2go(new_tax_id)
#match_uni_ids_from_list(uni_ids, newf)

ipr2go_f = "output/_bmq_ipr_go.txt" # Dictionary of ipr2go from EBI last release



def plot_enrichments():
    go_enr, ipr_enr, = enrichment_analysis_domains(choice, choice2, ipr2go_f)
    plot_2ylines2(ipr_enr)
    plot_2ylines2(go_enr)

####################################################################

def venn3(thefiles):
    files = ["theInputList1", "theInputList2", "theInputList3", "theInputList4","theInputList5","theInputList6","theInputList7","theOutput"]
    vars = []
    for key in files:
        #var = IntVar()
        vars.append(key)
    count = 1
    #print theFilenames1
    for i in range(len(thefiles[:3])):
        if len(thefiles[i]) > 1:
            vars[i] = thefiles[i]
            print i, vars[i]
        count += 1
    #theOutput = str(raw_input("Save file as : ")+"_merged.txt")
    get_domain_venn_diagram1(vars[0], vars[1], vars[2])

    #del theInputFiles[:]


###################################################################
### UPDATE FRAME
###################################################################


def togglex():
    if bt20.visible:
        btnTogglix["text"] = "Show "
        print "Now u don't"
        #mybutton.grid_remove()
        bt20.grid_remove()
        bt21.grid_remove()
        bt22.grid_remove()
        bt23.grid_remove()
        bt25.grid_remove()
        bt33.grid_remove()
        bt35.grid_remove()

        bt30.grid_remove()
        bt31.grid_remove()
        bt32.grid_remove()
        bt24.grid_remove()
        bt34.grid_remove()
        bt41.grid_remove()
        bt42.grid_remove()
        bt421.grid_remove()

        bt45.grid_remove()
        bt44.grid_remove()
        bt46.grid_remove()
        bt47.grid_remove()
        bt48.grid_remove()
        bt499.grid_remove()
        bt49.grid_remove()
        bt50.grid_remove()

        bt511.grid_remove()
        bt51.grid_remove() #bt11 = Button(frm1, text = " PLOT ENR ", command = plot_enrichments , bg="blue", fg="white").pack(padx=5, pady=5, side=LEFT)
        bt512.grid_remove()#bt11 = Button(frm1, text = " PLOT ENR ", command = plot_enrichments , bg="blue", fg="white").pack(padx=5, pady=5, side=LEFT)
        bt52.grid_remove()
        bt54.grid_remove()
        bt55.grid_remove()
        bt53.grid_remove()
        bt60.grid_remove()
        bt61.grid_remove()
        bt62.grid_remove()
        
    else:
        #mybutton.grid()
        bt20.grid(row=0, column = 0, sticky=N+S+E+W, padx=1, pady=1)
        bt21.grid(row=1, column = 0, sticky=N+S+E+W, padx=1, pady=1)
        bt22.grid(row=2, column = 0, sticky=N+S+E+W, padx=1, pady=1)
        bt23.grid(row=0, column = 1, sticky=N+S+E+W, padx=1, pady=1)
        bt25.grid(row=1, column = 1, sticky=N+S+E+W, padx=1, pady=1)
        bt33.grid(row=0, column = 9, sticky=N+S+E+W, padx=1, pady=1)
        bt35.grid(row=1, column = 9, sticky=N+S+E+W, padx=1, pady=1)

        bt30.grid(row=0, column = 2, sticky=N+S+E+W, padx=1, pady=1)
        bt31.grid(row=1, column = 2, sticky=N+S+E+W, padx=1, pady=1)
        bt32.grid(row=2, column = 2, sticky=N+S+E+W,padx=1, pady=1)
        bt24.grid(row=2, column = 1, sticky=N+S+E+W,padx=1, pady=1)
        bt34.grid(row=0, column = 3, sticky=N+S+E+W, padx=1, pady=1)
        bt41.grid(row=0, column = 4, sticky=N+S+E+W, padx=1, pady=1)
        bt42.grid(row=1, column = 3, sticky=N+S+E+W, padx=1, pady=1)
        bt421.grid(row=2, column = 3, sticky=N+S+E+W, padx=1, pady=1)


        bt45.grid(row=2, column = 4, sticky=N+S+E+W, padx=1, pady=1)
        bt44.grid(row=1, column = 4, sticky=N+S+E+W, padx=1, pady=1)
        bt46.grid(row=0, column = 5, sticky=N+S+E+W, padx=1, pady=1)
        bt47.grid(row=1, column = 5, sticky=N+S+E+W, padx=1, pady=1)
        bt48.grid(row=2, column = 5, sticky=N+S+E+W, padx=1, pady=1)
        bt499.grid(row=0, column = 6, sticky=N+S+E+W, padx=1, pady=1)
        bt49.grid(row=1, column = 6, sticky=N+S+E+W, padx=1, pady=1)
        bt50.grid(row=2, column = 6, sticky=N+S+E+W, padx=1, pady=1)

        bt511.grid(row=0, column = 7, sticky=N+S+E+W, padx=1, pady=1)
        bt51.grid(row=1, column = 7, sticky=N+S+E+W, padx=1, pady=1) #bt11 = Button(frm1, text = " PLOT ENR ", command = plot_enrichments , bg="blue", fg="white").pack(padx=5, pady=5, side=LEFT)
        bt512.grid(row=2, column = 7, sticky=N+S+E+W, padx=1, pady=1) #bt11 = Button(frm1, text = " PLOT ENR ", command = plot_enrichments , bg="blue", fg="white").pack(padx=5, pady=5, side=LEFT)
        bt52.grid(row=0, column = 8, sticky=N+S+E+W, padx=1, pady=1)
        bt54.grid(row=1, column = 8, sticky=N+S+E+W, padx=1, pady=1)
        bt55.grid(row=2, column = 8, sticky=N+S+E+W, padx=1, pady=1)
        bt53.grid(row=2, column = 9, sticky=N+S+E+W, padx=1, pady=1)
        bt60.grid(row=0, column = 10, sticky=N+S+E+W, padx=1, pady=1)
        bt61.grid(row=1, column = 10, sticky=N+S+E+W, padx=1, pady=1)
        bt62.grid(row=2, column = 10, sticky=N+S+E+W, padx=1, pady=1)

       

        print "Now usee it"
        btnTogglix["text"] = "Hide "
    #mybutton.visible = not mybutton.visible
    bt20.visible = not bt20.visible


#################
def togglixx():
    #def hide():
    if bt20.visible:
        bt20.withdraw()
        btnTogglix["text"] = "Show "
        print "Now u don't"
    
    else:

        
        #def show():
        btnTogglix["text"] = "Hide "
        print "u see it"
        bt20.update()
        bt20.deiconify()


    
# Update IMAGE

#img = ImageTk.PhotoImage(Image.open(path))
#panel = tk.Label(root, image = img)
#panel.pack(side = "bottom", fill = "both", expand = "yes")




###################################################################
################ FRAME 8
###################################################################



bt20= Button(frm8, text="LAcFilter", command= lambda: match_uni_ids_from_list(theUniIDs, theGeneIDs, choice) , bg="lightcyan2", fg="black") #bt20.pack( padx=5, pady=1, side=LEFT)
bt20.grid(row=0, column = 0, sticky=N+S+E+W, padx=1, pady=1)
bt20.visible = True

bt21 = Button(frm8, text="LIPRFilter", command= lambda:from_ipr_get_ac_infile(theIprIDs, choice) , bg="lightcyan2", fg="black") #bt21.pack( padx=5, pady=1, side=LEFT)
bt21.grid(row=1, column = 0, sticky=N+S+E+W, padx=1, pady=1)

bt22 = Button(frm8, text="LGOFilter", command= lambda:from_go_get_ac_infile(theGOIDs, choice) , bg="lightcyan2", fg="black") #bt22.pack( padx=5 , pady=1, side=LEFT)
bt22.grid(row=2, column = 0, sticky=N+S+E+W, padx=1, pady=1)

bt23 = Button(frm8, text="FindIn", command= lambda:find_word_in_line(choice,Words) , bg="lightcyan2", fg="black") #bt23.pack(padx=5, pady=1, side=LEFT)
bt23.grid(row=0, column = 1, sticky=N+S+E+W, padx=1, pady=1)

bt25 = Button(frm8, text ="Tax Name", command=lambda: from_name_get_tax_id(species_list), bg="lightcyan2" , fg="black") #bt25.pack(padx=5, pady=1, side=LEFT) # b1.pack(pady=1)
bt25.grid(row=1, column = 1, sticky=N+S+E+W, padx=1, pady=1)

bt33 = Button(frm8, text="TXT2EXCEL", command= lambda:convert_txt2excel(choice) , bg="red4", fg="white") #bt33.pack(padx=5, pady=1, side=LEFT)
bt33.grid(row=0, column = 9, sticky=N+S+E+W, padx=1, pady=1)

bt35 = Button(frm8, text="ExportAll", command= lambda:convertall_txt2excel(theInputFiles) , bg="red4", fg="white") #bt35.pack(padx=5, pady=1, side=LEFT)
bt35.grid(row=1, column = 9, sticky=N+S+E+W, padx=1, pady=1)

#bt550 = Button(frm8, text = " EnrichTacs ", command = lambda: enrichment_analysis_accessions(choice, choice2, ipr2go_f) , bg="blue", fg="white")
#bt550.pack(padx=5, pady=1,side=LEFT)


########################################################
################# FRAME 9    #   GRID BAR       ######## 
########################################################

def merge_me(choice, choice2):

    try:
        merge_ac_with_ipr_ipr(choice,choice2)
        #merge_inline(choice, choice2)
        
    except:
        pass
    


bt30 = Button(frm8, text="FAcFilter", command= lambda: match2files_by_ac(choice, choice2), bg="blue", fg="white") #bt30.pack(padx=5, pady=1, side=LEFT)
bt30.grid(row=0, column = 2, sticky=N+S+E+W, padx=1, pady=1)

bt31 = Button(frm8, text =" Merge by AC-IPR ", command=lambda:  merge_me(choice, choice2), bg="blue", fg="white") #bt31.pack(padx=5,  pady=1, side=LEFT) # check ok button, # save_unique(theInputFiles)
bt31.grid(row=1, column = 2, sticky=N+S+E+W, padx=1, pady=1)

    

bt32 = Button(frm8, text =" Match by AC", command=lambda: match2files_by_ac_only(choice, choice2) , bg="blue", fg="white") #bt32.pack(padx=5, pady=1, side=LEFT) # check ok button, # save_unique(theInputFiles)
bt32.grid(row=2, column = 2, sticky=N+S+E+W,padx=1, pady=1)


bt24 = Button(frm8, text=" Join ", command= lambda: join_lines(choice,choice2) , bg="blue", fg="white") #bt24.pack(padx=5, pady=1, side=LEFT)
bt24.grid(row=2, column = 1, sticky=N+S+E+W,padx=1, pady=1)

#child2family(sourceFile, filterFile)

bt34 = Button(frm8, text=" SUM ", command= lambda:readUniProtTab_sum(choice, n = 0, with_fields = False) , bg="blue", fg="white") #bt34.pack(padx=5, pady=1, side=LEFT)
bt34.grid(row=0, column = 3, sticky=N+S+E+W, padx=1, pady=1)



bt41 = Button(frm8, text ="Bar Plot", command=lambda: bar_plot2labels(theInputFiles), bg="DeepSkyBlue3", fg="white")  #bt41.pack(padx=5,  pady=1, side=LEFT)# ()print_data
bt41.grid(row=0, column = 4, sticky=N+S+E+W, padx=1, pady=1)

bt42 =Button(frm8, text ="Scatter GO", command=lambda: scatter_go_functions(choice), bg="DeepSkyBlue3", fg="white") #bt42.pack(padx=5, pady=1,side=LEFT)# ()print_data
bt42.grid(row=1, column = 3, sticky=N+S+E+W, padx=1, pady=1)


bt421 =Button(frm8, text ="Scatter IPR", command=lambda: scatter_ipr_functions(choice), bg="DeepSkyBlue3", fg="white") #bt42.pack(padx=5, pady=1,side=LEFT)# ()print_data
bt421.grid(row=2, column = 3, sticky=N+S+E+W, padx=1, pady=1)


#bt422 =Button(frm8, text ="Scatter PATH", command=lambda: scatter_path(choice), bg="blue", fg="white") #bt42.pack(padx=5, pady=1,side=LEFT)# ()print_data
#bt422.grid(row=1, column = 3, sticky=N+S+E+W, padx=1, pady=1)

#bb22 = Button(frm1, text ="Summary Stats", command= lambda: summary_funtions_1(theFilenames)).pack(side=LEFT) # check ok button
#b2 = Button(frm1, text ="Summary Analysis", command=summary_analysis).pack(side=LEFT) # check ok button

#"from auxiliary_functions_tk import *
# "IPR2GO_bsbmv" = lambda: from_ac2iprbs_get_ipr2gobmv("5671_bmq_ac_ipr.txt", "bmv_ipr_go.txt", "5671_bmq_ac_ipr_go44.txt")) # check ok button


#Button(frm1,text='Enter Tax ID', command=ask_tax_id).pack(fill=X)


######################################################################################

def validate_venn2():

    global choice
    global choice2
    if choice in theInputFiles and choice2 in theInputFiles:
        read_tab_enrichment_ipr(choice, choice2)
    else:
        showinfo("Warning!", "Combobox (2) not valid file!")

import run_biopython_tk
import run_bioblast_tk


bt44 = Button(frm8, text ="Pie IPR", command=lambda : count_ipr_pie(choice, ipr2go_f) , bg="DeepSkyBlue3", fg="white") #bt44.pack(padx=5, pady=1, side=LEFT)# ()print_data
bt44.grid(row=1, column = 4, sticky=N+S+E+W, padx=1, pady=1)

bt45 = Button(frm8, text ="Pie GO", command=lambda : count_go_pie(choice, ipr2go_f) , bg="DeepSkyBlue3", fg="white") #bt45.pack(padx=5,  pady=1,side=LEFT)# ()print_data
bt45.grid(row=2, column = 4, sticky=N+S+E+W, padx=1, pady=1)

bt46 = Button(frm8, text ="GO CC", command= lambda: getGOcomponent(choice), bg="DeepSkyBlue3", fg="white") #bt46.pack(padx=5, pady=1, side=LEFT) # check ok button
bt46.grid(row=0, column = 5, sticky=N+S+E+W, padx=1, pady=1)


bt47 = Button(frm8, text ="GO BP", command= lambda: getGOprocess(choice), bg="DeepSkyBlue3", fg="white") # , cursor="hand2" # bt47.pack(padx=5, pady=1, side=LEFT) # check ok button
bt47.grid(row=1, column = 5, sticky=N+S+E+W, padx=1, pady=1)

bt48 = Button(frm8, text ="GO MF", command= lambda: getGOfunction(choice), bg="DeepSkyBlue3", fg="white") # bt48.pack(padx=5, pady=1, side=LEFT) # check ok button
bt48.grid(row=2, column = 5, sticky=N+S+E+W, padx=1, pady=1)


bt499 = Button(frm8, text ="AC Venn2", command=lambda: read_tab_enrichment_acs(choice, choice2), bg="DeepSkyBlue3", fg="white") #bt49.pack(padx=5, pady=1,side=LEFT) # check ok button, venn_diagram_manual
bt499.grid(row=0, column = 6, sticky=N+S+E+W, padx=1, pady=1)


bt49 = Button(frm8, text ="IPR Venn2", command=lambda: read_tab_enrichment_ipr(choice, choice2), bg="DeepSkyBlue3", fg="white") #bt49.pack(padx=5, pady=1,side=LEFT) # check ok button, venn_diagram_manual
bt49.grid(row=1, column = 6, sticky=N+S+E+W, padx=1, pady=1)


bt50 = Button(frm8, text ="IPR Venn3", command=lambda: venn3(theInputFiles) , bg="DeepSkyBlue3", fg="white") #bt50.pack(padx=5, pady=1,side=LEFT) # check ok button, venn_diagram_manual
bt50.grid(row=2, column = 6, sticky=N+S+E+W, padx=1, pady=1) # lambda: venn3(theInputFiles);pop_up_venn3_files




#####################################################

#b55 = Button(frm9, text ="BioBlast", command= lambda: access_bioblast_with_biopython1(theUniIDs), bg="blue", fg="white").pack(padx=5, pady=5, side=LEFT) # check ok button:access_blast_with_biopython,access_bioblast_with_biopython1

#b55 = Button(frm9, text ="BioBlast", command= lambda: access_blast_with_biopython3(theUniIDs), bg="blue", fg="white").pack(padx=5, pady=5, side=LEFT) # check ok button:access_blast_with_biopython,access_bioblast_with_biopython1

# BioBlast Homologues . ImportError: DLL load failed: %1 is not a valid Win32 application.

#####################################################


bt511 = Button(frm8, text = " EnrichPROT ", command = lambda: enrichment_analysis_proteins(choice, choice2, ipr2go_f) , bg="green3", fg="white") #bt511.pack(padx=5, pady=1,side=LEFT)
bt511.grid(row=0, column = 7, sticky=N+S+E+W, padx=1, pady=1)

bt51 = Button(frm8, text = " EnrichIPRGO ", command = lambda: enrichment_analysis_domains(choice, choice2, ipr2go_f) , bg="green3", fg="white") #bt51.pack(padx=5, pady=1,side=LEFT)
bt51.grid(row=1, column = 7, sticky=N+S+E+W, padx=1, pady=1) #bt11 = Button(frm1, text = " PLOT ENR ", command = plot_enrichments , bg="blue", fg="white").pack(padx=5, pady=5, side=LEFT)


bt512 = Button(frm8, text = " EnrichPATH ", command = lambda: enrichment_analysis_paths(choice, choice2, ipr2go_f) , bg="green3", fg="white") #bt51.pack(padx=5, pady=1,side=LEFT)
bt512.grid(row=2, column = 7, sticky=N+S+E+W, padx=1, pady=1) #bt11 = Button(frm1, text = " PLOT ENR ", command = plot_enrichments , bg="blue", fg="white").pack(padx=5, pady=5, side=LEFT)

bt52 = Button(frm8, text = "Enrich Pvalue", command = lambda: plot_2ylines_pval(choice) , bg="green4", fg="white") #bt52.pack(padx=5, pady=1, side=LEFT)
bt52.grid(row=0, column = 8, sticky=N+S+E+W, padx=1, pady=1)

bt54 = Button(frm8, text = "Enrich OR", command = lambda: plot_2ylines_or(choice) , bg="green4", fg="white") #bt54.pack(padx=5, pady=1, side=LEFT)
bt54.grid(row=1, column = 8, sticky=N+S+E+W, padx=1, pady=1)

bt55 = Button(frm8, text = " FDR 35", command = lambda: plot_2ylines_fdr(choice) , bg="green4", fg="white") #bt55.pack(padx=5, pady=1, side=LEFT)
bt55.grid(row=2, column = 8, sticky=N+S+E+W, padx=1, pady=1)

bt53 = Button(frm8, text = "FDR FWER", command = lambda: plot_2ylines_corr(choice) , bg="green4", fg="white") #bt53.pack(padx=5, pady=1, side=LEFT)
bt53.grid(row=2, column = 9, sticky=N+S+E+W, padx=1, pady=1)

bt60 = Button(frm8, text = "BioBlast", command = lambda: access_bioblast_with_biopython1(theUniIDs) , bg="lightcyan2", fg="black") #bt53.pack(padx=5, pady=1, side=LEFT)
bt60.grid(row=0, column = 10, sticky=N+S+E+W, padx=1, pady=1)

bt61 = Button(frm8, text = "GO Matrix", command = lambda: access_bioblast_go_matrix(theUniIDs) , bg="lightcyan2", fg="black") #bt53.pack(padx=5, pady=1, side=LEFT)
bt61.grid(row=1, column = 10, sticky=N+S+E+W, padx=1, pady=1)

bt62 = Button(frm8, text = "Abstracts", command = lambda: access_bioblast_abstracts(theUniIDs) , bg="lightcyan2", fg="black") #bt53.pack(padx=5, pady=1, side=LEFT)
bt62.grid(row=2, column = 10, sticky=N+S+E+W, padx=1, pady=1)




######################################################################
#####     STATUS BAR
#######################################################################33

messag20 = " Use ListBox to Filter AC in Selected File. "
messag21 = " Use ListBox to Filter IPR in Selected File. "
messag22 = " Use ListBox to Filter GO in Selected File. "
messag23 = " Use ListBox to Find Terms in Selected File. "
messag24 = " Use Comboboxes to join 2 files. "
messag25 = " Use ListBox to Enter a list of Tax Names. "

messag30 = " Use ComboBoxes to Filter by AC. "
messag31 = " Use ComboBoxes to Merge, same line, by AC. "
messag32 = " Use ComboBoxes to Match, same column, by AC. "
messag33 = " Use ComboBox (1) to Convert .txt file to .xls file."
messag34 = " Use ComboBox (1) to Count all different hits in column. "
messag35 = " Export all files from ComboBox, in EXCEL sheets."
messag36 = ""


#messag40 = "Use ComboBox to "
messag41 = " Use ComboBox (1) to Plot Bars ( all files ) of AC-IPR-GO. "
messag42 = " Use ComboBox (1) to Plot Scatter of GO Annotations. "

#messag43 = ""
messag44 = " Use ComboBox (1) to Plot Pie of IPR Domain Functions. "
messag45 = " Use ComboBox (1) to Plot Pie of GO Annotations. "

messag46 = " Use ComboBox (1) to Plot Pie of GO Cellular Component. Bio:SwissProt"
messag47 = " Use ComboBox (1) to Plot Pie of GO Biological Process. Bio:SwissProt"
messag48 = " Use ComboBox (1) to Plot Pie of GO Molecular Function. Bio:SwissProt"

messag499 = "Use ComboBoxes to Plot Venn2 of ACs."
messag49 = " Use ComboBoxes to Plot Venn2 of IPR Domain Functions. "
messag50 = " Use ComboBox (1) to Plot Venn3 ( first 3 files ) of IPR Domains.  "

messag511 = " Use ComboBoxes to Calculate Enrichment of IPR / GO  annotations in ACsIPR. "
messag510 = " Use ComboBoxes to Calculate Enrichment of IPR / GO  annotations in total IDs. "
messag512 = " Use ComboBoxes to Calculate Enrichment of IPR / GO  annotations in unique IDs. "
messag52 = " Use ComboBox (1) to Plot of IPR / GO Functions versus P-Value. "
messag53 = " Use ComboBox (1) to Plot Graph of IPR / GO Domain Functions versus FDR | FWER. "
messag54 = " Use ComboBox (1) to Plot Graph of IPR / GO Domain Functions versus Odds Ratio. "

messag55 = " Use ComboBox (1) to Plot Graph of IPR / GO Domain Functions versus FDR. "


messag60 = " Use ComboBox (1) to retrieve protein homologues . "
messag61 = " Use ComboBox (1) to construct GO Matrix . "
messag62 = " Use ComboBox (1) to retrieve Abstracts . "


########################
def on_enter( event):
    l2.configure(text="Hello World!")


def on_leave( enter):
    l2.configure(text="")

#########################

def on_enter20( event):
    l2.configure(text=messag20)

def on_enter21( event):
    l2.configure(text=messag21)

def on_enter22( event):
    l2.configure(text=messag22)

def on_enter23( event):
    l2.configure(text=messag23)


def on_enter24( event):
    l2.configure(text=messag24)

def on_enter25( event):
    l2.configure(text=messag25)
#########################

def on_enter30( event):
    l2.configure(text=messag30)

def on_enter31( event):
    l2.configure(text=messag31)

def on_enter32( event):
    l2.configure(text=messag32)

def on_enter33( event):
    l2.configure(text=messag33)

def on_enter34( event):
    l2.configure(text=messag34)

def on_enter35( event):
    l2.configure(text=messag35)
#########################

def on_enter40( event):
    l2.configure(text=messag40)

def on_enter41( event):
    l2.configure(text=messag41)

def on_enter42( event):
    l2.configure(text=messag42)

def on_enter43( event):
    l2.configure(text=messag43)

def on_enter44( event):
    l2.configure(text=messag44)

def on_enter45( event):
    l2.configure(text=messag45)

def on_enter46( event):
    l2.configure(text=messag46)

def on_enter47( event):
    l2.configure(text=messag47)

def on_enter48( event):
    l2.configure(text=messag48)

def on_enter49( event):
    l2.configure(text=messag49)

def on_enter499(event):
    l2.configure(text=messag499)

def on_enter50( event):
    l2.configure(text=messag50)

def on_enter511( event):
    l2.configure(text=messag511)


def on_enter510( event):
    l2.configure(text=messag510)

def on_enter512( event):
    l2.configure(text=messag512)
    
def on_enter52( event):
    l2.configure(text=messag52)

def on_enter53( event):
    l2.configure(text=messag53)

def on_enter54( event):
    l2.configure(text=messag54)

def on_enter55( event):
    l2.configure(text=messag55)
    

def on_enter60( event):
    l2.configure(text=messag60)

def on_enter61( event):
    l2.configure(text=messag61)

def on_enter62( event):
    l2.configure(text=messag62)


############################  STATUS BAR  ########################
frm17 = Frame(jp)

frm17.pack(side="top", fill="both", expand="true", pady=3)

# MAKE SHOW/HIDE BUTTON FOR GRID BUTTONS
btnTogglix = Button(frm17, text="Hide ", command=togglex)
btnTogglix.pack(side=LEFT)


l1 = Label(frm17, text=" Status ")# hover over me
l1.pack(side=LEFT)

l2 = Label(frm17,text="", width=100)
l2.pack(side=LEFT, fill="x")


l1.bind("<Enter>", on_enter)
l1.bind("<Leave>", on_leave)


###############################

bt20.bind("<Enter>", on_enter20)
bt20.bind("<Leave>", on_leave)

bt21.bind("<Enter>", on_enter21)
bt21.bind("<Leave>", on_leave)

bt22.bind("<Enter>", on_enter22)
bt22.bind("<Leave>", on_leave)

bt23.bind("<Enter>", on_enter23)
bt23.bind("<Leave>", on_leave)


bt24.bind("<Enter>", on_enter24)
bt24.bind("<Leave>", on_leave)

bt25.bind("<Enter>", on_enter25)
bt25.bind("<Leave>", on_leave)

###############################

bt30.bind("<Enter>", on_enter30)
bt30.bind("<Leave>", on_leave)

bt31.bind("<Enter>", on_enter31)
bt31.bind("<Leave>", on_leave)

bt32.bind("<Enter>", on_enter32)
bt32.bind("<Leave>", on_leave)

bt33.bind("<Enter>", on_enter33)
bt33.bind("<Leave>", on_leave)

bt34.bind("<Enter>", on_enter34)
bt34.bind("<Leave>", on_leave)

bt35.bind("<Enter>", on_enter35)
bt35.bind("<Leave>", on_leave)

###############################

#bt40.bind("<Enter>", on_enter40)
#bt40.bind("<Leave>", on_leave)

bt41.bind("<Enter>", on_enter41)
bt41.bind("<Leave>", on_leave)

bt42.bind("<Enter>", on_enter42)
bt42.bind("<Leave>", on_leave)

#bt43.bind("<Enter>", on_enter43)
#bt43.bind("<Leave>", on_leave)

bt44.bind("<Enter>", on_enter44)
bt44.bind("<Leave>", on_leave)

bt45.bind("<Enter>", on_enter45)
bt45.bind("<Leave>", on_leave)

bt46.bind("<Enter>", on_enter46)
bt46.bind("<Leave>", on_leave)

bt47.bind("<Enter>", on_enter47)
bt47.bind("<Leave>", on_leave)

bt48.bind("<Enter>", on_enter48)
bt48.bind("<Leave>", on_leave)

bt49.bind("<Enter>", on_enter49)
bt49.bind("<Leave>", on_leave)

bt499.bind("<Enter>", on_enter499)
bt499.bind("<Leave>", on_leave)

bt50.bind("<Enter>", on_enter50)
bt50.bind("<Leave>", on_leave)

bt511.bind("<Enter>", on_enter511)
bt511.bind("<Leave>", on_leave)

bt51.bind("<Enter>", on_enter510)
bt51.bind("<Leave>", on_leave)

bt52.bind("<Enter>", on_enter512)
bt52.bind("<Leave>", on_leave)

bt52.bind("<Enter>", on_enter52)
bt52.bind("<Leave>", on_leave)

bt53.bind("<Enter>", on_enter53)
bt53.bind("<Leave>", on_leave)

bt54.bind("<Enter>", on_enter54)
bt54.bind("<Leave>", on_leave)

bt55.bind("<Enter>", on_enter55)
bt55.bind("<Leave>", on_leave)


bt60.bind("<Enter>", on_enter60)
bt60.bind("<Leave>", on_leave)

bt61.bind("<Enter>", on_enter61)
bt61.bind("<Leave>", on_leave)

bt62.bind("<Enter>", on_enter62)
bt62.bind("<Leave>", on_leave)
#################################






#statbar(frm17).pack(side="top", fill="both", expand="true")

#Example(frm17).pack(side="top", fill="both", expand="true")




######################################################################

import winsound
# ValueError: frequency must be in 37 thru 32767
Freq = 2500 # Set Frequency To 2500 Hertz
Dur = 750 # Set Duration To 1000 ms == 1 second
winsound.Beep(Freq,Dur)

###################################################################3
def update_timeText():
    if (state):
        global timer
        # Every time this function is called,
        # we will increment 1 centisecond (1/100 of a second)
        timer[2] += 1
        # Every 100 centisecond is equal to 1 second
        if (timer[2] >= 100):
            timer[2] = 0
            timer[1] += 1
        # Every 60 seconds is equal to 1 min
        if (timer[1] >= 60):
            timer[0] += 1
            timer[1] = 0
        # We create our time string here
        timeString = pattern.format(timer[0], timer[1], timer[2])
        # Update the timeText Label box with the current time
        timeText.configure(text=timeString)
        # Call the update_timeText() function after 1 centisecond
    jp.after(10, update_timeText)





###################################################################3
#######      TIME  STATUS
##################################################3

# Simple status flag
# False mean the timer is not running
# True means the timer is running (counting)
state = True

# Our time structure [min, sec, centsec]
timer = [0, 0, 0]
# The format is padding all the
pattern = '{0:02d}:{1:02d}:{2:02d}'

# Create a timeText Label (a text box)

timeText = Label(frm17, text=" 00:00:00 ", font=("Helvetica", 10))
timeText.pack(side = LEFT)
update_timeText()



#####################################################################

"""with open('mydata.txt') as fp:
    for line in iter(fp.readline, ''):
        process_line(line)

def ensure_dir(f):
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)
        """

##################################################################333



time_end = datetime.datetime.fromtimestamp(time.time())
print("Time elapsed: ", str(time_end - time_begin))


frm20 = Frame(jp)
frm20.pack(side="top", fill="both", expand="true")
Button(frm20, text='Quit', command=get_out, fg="blue", font=('arial',11, 'bold')).pack(fill=X , expand = "True", side=TOP)


jp.mainloop()


######################################################################3
## Pathosystems - *ATRIC
## http://patricbrc.org/portal/portal/patric/GenomeList?cType=taxon&cId=2&dataSource=&displayMode=&pk=&kw=

## NCBI
## http://www.ncbi.nlm.nih.gov/protein/?term=txid435258[Organism:noexp]
# jpcm5
## http://www.ncbi.nlm.nih.gov/protein/?term=txid5671%5BOrganism%3Anoexp%5D
# Leish infantum

## http://www.ncbi.nlm.nih.gov/gene/?term=txid435258[Organism:noexp]
# total no.8381
## http://www.ncbi.nlm.nih.gov/proteinclusters/?term=txid435258[Organism:noexp]
# total no.7008
# Bioython cookbook
# http://biopython.org/DIST/docs/tutorial/Tutorial.html

# OPEN SANGER SPECIE
# http://www.sanger.ac.uk/resources/downloads/protozoa/leishmania-infantum.html

# Genome information by organism
# http://www.ncbi.nlm.nih.gov/genome/browse/

# quick guide genomes sequenced
# http://www.genomenewsnetwork.org/resources/sequenced_genomes/genome_guide_p1.shtml

# 26570 genomes found
# http://patricbrc.org/portal/portal/patric/GenomeList?cType=taxon&cId=2

# National Microbial Pathogens Data Resource
# http://www.nmpdr.org/FIG/wiki/view.cgi/FIG/TaxonomyIdentifier

# NCBI Browser
# http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Root

# HELP List
# Ensembl REST API Endpoints
# http://rest.ensembl.org/

# BROWSE TAX ID
# http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/

# http://www.genome.jp/kegg/

# http://www.programmableweb.com/api/ebi-interproscan
# http://www.ebi.ac.uk/interpro/about.html

# https://www.python.org/download/releases/2.7.3/
## https://pythonhosted.org/bioservices/
## http://www.ebi.ac.uk/Tools/webservices/services/archive/pfa/iprscan_soap
## http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_rest
## http://www.ebi.ac.uk/Tools/webservices/services/pfa/iprscan5_soap
## http://www.ebi.ac.uk/Tools/webservices/ = FULL LIST OF JOBS
#

# LESMA E O GENE DA LESMA MARINHA
# Elysia Chlorotica feeds on Vaucheria litorea
# http://www.jornalciencia.com/meio-ambiente/animais/4638-animal-e-vegetal-ao-mesmo-tempo-lesma-marinha-faz-fotossintese-apos-roubar-genes-de-algas
#
#######################################################################
# Print, # SHOW # Make  # Quit # Clear # Run File, save, submit

# Status # Run Code ###  Legenda # GRAPHICS


# Access database : http://www.tutorialspoint.com/python/python_database_access.htm

# tax ids:http://pythonhosted.org/Orange-Bioinformatics/reference/taxonomy.html

# openurl: http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/
# http://www.uniprot.org/help/taxonomy
#


# ftp://ftp.ebi.ac.uk/pub/databases/interpro/names.dat
# ftp://ftp.ebi.ac.uk/pub/databases/interpro/entry.list
