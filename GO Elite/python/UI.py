###UI
#Copyright 2005-2008 J. David Gladstone Institutes, San Francisco California
#Author Nathan Salomonis - nsalomonis@gmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

"""This module contains fairly generic routines for reading in user interface (UI) instructions from
an existing configuration options file (options.txt), buiding those interfaces, downloading and/or
processing files."""

import math
import statistics
import sys, string
import os.path, platform
import time
import webbrowser
import shutil
import update; reload(update)
import BuildAffymetrixAssociations; reload(BuildAffymetrixAssociations)
import gene_associations; reload(gene_associations)
import export
import unique
import OBO_import
import datetime
import traceback

try:
    import WikiPathways_webservice
except Exception:
    print 'WikiPathways visualization not supported (requires installation of suds)'
try:
    from PIL import Image as PIL_Image
    import ImageTk
except Exception:
    print 'Python Imaging Library not installed... using default PNG viewer'
from sys import argv

try:
    import Tkinter
    from Tkinter import *
    import PmwFreeze
    from Tkconstants import LEFT
    import tkMessageBox
    import tkFileDialog
except ImportError: print "\nPmw or Tkinter not found... proceeding with manual input"

mac_print_mode = 'no' 
if os.name == 'posix':  mac_print_mode = 'yes' #os.name is  'posix', 'nt', 'os2', 'mac', 'ce' or 'riscos'
debug_mode = 'no'

def filepath(filename):
    fn = unique.filepath(filename)
    return fn

def osfilepath(filename):
    fn = filepath(filename)
    fn = string.replace(fn,'\\','/')
    return fn

def read_directory(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt" or entry[-4:] == ".csv" or ".zip" in entry or '.obo' in entry or '.ontology' in entry: dir_list2.append(entry)
    return dir_list2

def readDirText(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Code to prevent folder names from being included
    for entry in dir_list:
        if entry[-4:] == ".txt": dir_list2.append(entry)
    return dir_list2

def getFolders(sub_dir):
    dir_list = unique.read_directory(sub_dir); dir_list2 = []
    ###Only get folder names
    for entry in dir_list:
        if entry[-4:] != ".txt" and entry[-4:] != ".csv" and ".zip" not in entry: dir_list2.append(entry)
    return dir_list2

def returnDirectoriesNoReplace(dir):
    dir_list = unique.returnDirectoriesNoReplace(dir); dir_list2 = []
    for entry in dir_list:
        if '.' not in entry: dir_list2.append(entry)
    return dir_list2

def returnFilesNoReplace(dir):
    dir_list = unique.returnDirectoriesNoReplace(dir); dir_list2 = []
    for entry in dir_list:
        if '.' in entry: dir_list2.append(entry)
    return dir_list2
    
def cleanUpLine(line):
    line = string.replace(line,'\n','')
    line = string.replace(line,'\c','')
    data = string.replace(line,'\r','')
    data = string.replace(data,'"','')
    return data

################# GUI #################

def wait(seconds):
    ### stalls the analysis for the designated number of seconds (allows time to print something to the GUI)
    start_time = time.time()
    diff = 0
    while diff<seconds:
        diff = time.time()-start_time
        
class GUI:
    def ViewWikiPathways(self):
        """ Canvas is already drawn at this point from __init__ """
        global pathway_db
        pathway_db={}
        button_text = 'Help'

        ### Create a species drop-down option that can be updated
        current_species_names = getSpeciesList()
        self.title = 'Select species to search for WikiPathways '
        self.option = 'species_wp'
        self.options = ['---']+current_species_names #species_list
        self.default_option = 0
        self.dropDown()
        
        ### Create a label that can be updated below the dropdown menu
        self.label_name = StringVar()
        self.label_name.set('Pathway species list may take several seconds to load')
        self.invokeLabel() ### Invoke a new label indicating that the database is loading
                    
        ### Create a MOD selection drop-down list
        null,system_list,mod_list = importSystemInfo()
        self.title = 'Select the ID system to translate to (MOD)'
        self.option = 'mod_wp'
        self.options = mod_list
        try: self.default_option = mod_list.index('Ensembl') ### Get the Ensembl index number
        except Exception: self.default_option = 0
        self.dropDown()
        
        ### Create a file selection option
        self.title = 'Select GO-Elite input ID text file'
        self.notes = 'note: ID file must have a header row and at least three columns:\n'
        self.notes += '(1) Identifier, (2) System Code, (3) Value to map (- OR +)\n'
        self.file_option = 'goelite_input_file'
        self.directory_type = 'file'
        self.FileSelectionMenu()
        
        dispaly_pathway = Button(text = 'Display Pathway', command = self.displayPathway)
        dispaly_pathway.pack(side = 'right', padx = 10, pady = 10)

        back_button = Button(self._parent, text="Back", command=self.goBack) 
        back_button.pack(side = 'right', padx =10, pady = 5)
        
        quit_win = Button(self._parent, text="Quit", command=self.quit) 
        quit_win.pack(side = 'right', padx =10, pady = 5)

        try: help_button = Button(self._parent, text=button_text, command=self.GetHelpTopLevel); help_button.pack(side = 'left', padx = 5, pady = 5)
        except Exception: help_button = Button(self._parent, text=button_text, command=self.linkout); help_button.pack(side = 'left', padx = 5, pady = 5)

        self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        self._parent.mainloop()

    def FileSelectionMenu(self):
        option = self.file_option
        group = PmwFreeze.Group(self.parent_type,tag_text = self.title)
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 2)
        
        def filecallback(callback=self.callback,option=option): self.getPath(option)
        default_option=''
        entrytxt = StringVar(); #self.entrytxt.set(self.default_dir)
        entrytxt.set(default_option)
        self.pathdb[option] = entrytxt
        self._user_variables[option] = default_option
        entry = Entry(group.interior(),textvariable=self.pathdb[option]);
        entry.pack(side='left',fill = 'both', expand = 0.7, padx = 10, pady = 2)
        button = Button(group.interior(), text="select "+self.directory_type, width = 10, fg="red", command=filecallback)
        button.pack(side=LEFT, padx = 2,pady = 2)
        if len(self.notes)>0: ln = Label(self.parent_type, text=self.notes,fg="blue"); ln.pack(padx = 10)
        
    def dropDown(self):
        def comp_callback(tag,callback=self.callbackWP,option=self.option):
            callback(tag,option)
        self.comp = PmwFreeze.OptionMenu(self.parent_type,
            labelpos = 'w', label_text = self.title, items = self.options, command = comp_callback)
        if self.option == 'wp_id_selection':
            self.wp_dropdown = self.comp ### update this variable later (optional)
        self.comp.pack(anchor = 'w', padx = 10, pady = 0)
        self.comp.invoke(self.default_option) ###Just pick the first option

    def comboBox(self):
        """ Alternative, more sophisticated UI than dropDown (OptionMenu).
        Although it behaves similiar it requires different parameters, can not be
        as easily updated with new lists (different method) and requires explict
        invokation of callback when a default is set rather than selected. """
        
        def comp_callback(tag,callback=self.callbackWP,option=self.option):
            callback(tag,option)
        self.comp = PmwFreeze.ComboBox(self.parent_type,
            labelpos = 'w', dropdown=1, label_text = self.title,
            unique = 0, history = 0,
            scrolledlist_items = self.options, selectioncommand = comp_callback)
        
        try: self.comp.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.comp.invoke())
        except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
        
        if self.option == 'wp_id_selection':
            self.wp_dropdown = self.comp ### update this variable later (optional)
        self.comp.pack(anchor = 'w', padx = 10, pady = 0)
        self.comp.selectitem(self.default_option) ###Just pick the first option
        self.callbackWP(self.options[0],self.option)  ### Explicitly, invoke first option (not automatic)
        
    def invokeLabel(self):
        self.label_object = Label(self.parent_type, textvariable=self.label_name,fg="blue"); self.label_object.pack(padx = 10)
        
    def invokeStatusLabel(self):
        self.label_object = Label(self.parent_type, textvariable=self.label_status_name,fg="blue"); self.label_object.pack(padx = 10)
          
    def enterMenu(self):
        if len(self.notes)>0:
            lb = Label(self.parent_type, text=self.notes,fg="black"); lb.pack(pady = 5)
        ### Create and pack a horizontal RadioSelect widget
        def custom_validate(tag,custom_validate=self.custom_validate,option=self.option):
            validate = custom_validate(tag,self.option)
        self.entry_field = PmwFreeze.EntryField(self.parent_type,
                labelpos = 'w', label_text = self.title, validate = custom_validate, 
                value = self.default_option, hull_borderwidth = 2)
        self.entry_field.pack(fill = 'x', expand = 0.7, padx = 10, pady = 5)
    
    def displayPathway(self):
        filename = self._user_variables['goelite_input_file']
        mod_type = self._user_variables['mod_wp']
        species = self._user_variables['species_wp']
        pathway_name = self._user_variables['wp_id_selection']
        wpid_selected = self._user_variables['wp_id_enter']
        species_code = species_codes[species].SpeciesCode()
        wpid = None
        if len(wpid_selected)>0:
            wpid = wpid_selected
        elif len(self.pathway_db)>0:
            for wpid in self.pathway_db:
                if pathway_name == self.pathway_db[wpid].WPName():
                    break
        if len(filename)==0:
            print_out = 'Select an input ID file with values first'
            WarningWindow(print_out,'Error Encountered!')
        else:
            try:
                self.graphic_link = WikiPathways_webservice.visualizePathwayAssociations(filename,species_code,mod_type,wpid)
                self.wp_status = 'Pathway images colored and saved to disk by webservice\n(see image title for location)'
                self.label_status_name.set(self.wp_status)
                try: self.viewPNGFile() ### ImageTK PNG viewer
                except Exception:
                    try: self.openPNGImage() ### OS default PNG viewer
                    except Exception:
                        self.wp_status = 'Unable to open PNG file using operating system'
                        self.label_status_name.set(self.wp_status)
            except Exception,e:
                try:
                    wp_logfile = filepath('webservice.log')
                    wp_report = open(wp_logfile,'a')
                    wp_report.write(traceback.format_exc())
                except Exception:
                    None
                try:
                    print traceback.format_exc()
                except Exception:
                    null=None ### Occurs when transitioning back from the Official Database download window (not sure why) -- should be fixed in 1.2.5 (sys.stdout not re-routed)
                if 'force_no_matching_error' in traceback.format_exc():
                    print_out = 'None of the input IDs mapped to this pathway'
                elif 'force_invalid_pathway' in traceback.format_exc():
                    print_out = 'Invalid pathway selected'
                elif 'IndexError' in traceback.format_exc():
                    print_out = 'Input ID file does not have at least 3 columns, with the second column being system code'
                elif 'ValueError' in traceback.format_exc():
                    print_out = 'Input ID file error. Please check that you do not have extra rows with no data' 
                elif 'source_data' in traceback.format_exc():
                    print_out = 'Input ID file does not contain a valid system code' 
                else:
                    print_out = 'Error generating the pathway "%s"' % pathway_name
                WarningWindow(print_out,'Error Encountered!')
        
    def getSpeciesPathways(self,species_full):
        pathway_list=[]
        self.pathway_db = WikiPathways_webservice.getAllSpeciesPathways(species_full)
        for wpid in self.pathway_db:
            pathway_list.append(self.pathway_db[wpid].WPName())
        pathway_list = unique.unique(pathway_list)
        pathway_list.sort()
        return pathway_list

    def callbackWP(self, tag, option):
        #print 'Button',[option], tag,'was pressed.'
        self._user_variables[option] = tag
        if option == 'species_wp':
            ### Add additional menu options based on user selection
            if tag != '---':
                ### If this already exists from an earlier iteration
                hault = False
                self.label_name.set('Loading available WikiPathways')
                try:
                    self.pathway_list=self.getSpeciesPathways(tag)
                    traceback_printout = ''
                except Exception,e:
                    if 'not supported' in traceback.format_exc():
                        print_out = 'Species not available at WikiPathways'
                        WarningWindow(print_out,'Species Not Found!')
                        traceback_printout=''
                        hault = True
                    elif 'URLError' in traceback.format_exc():
                        print_out = 'Internet connection could not be established'
                        WarningWindow(print_out,'Internet Error')
                        traceback_printout=''
                        hault = True
                    else:
                        traceback_printout = traceback.format_exc()
                    try: 
                        if len(self.pathway_list)>0: ### When true, a valid species was selected in a prior interation invoking the WP fields (need to repopulate)
                            hault = False
                    except Exception: None
                    self.pathway_list = ['None']; self.pathway_db={}
                self.label_name.set('')
                if hault == False:
                    try:
                        ### If the species specific wikipathways drop down exists, just update it
                        self.wp_dropdown._list.setlist(self.pathway_list)
                        self.wp_dropdown.selectitem(self.pathway_list[0])
                        self.callbackWP(self.pathway_list[0],'wp_id_selection')
                    except Exception:
                        ### Create a species specific wikipathways drop down
                        self.option = 'wp_id_selection'
                        self.title = 'Select WikiPathways to visualize your data'
                        if len(traceback_printout)>0:
                            self.title += traceback_printout ### Display the actual problem in the GUI (sloppy but efficient way for users to indicate the missing driver)
                        self.options = self.pathway_list
                        self.default_option = 0
                        self.comboBox() ### Better UI for longer lists of items (dropDown can't scroll on Linux)
                        
                        ### Create a species specific wikipathways ID enter option
                        self.notes = 'OR'
                        self.option = 'wp_id_enter'
                        self.title = 'Enter the WPID (example: WP254) '
                        self.default_option = ''
                        self.enterMenu()
                    try:
                        ### Create a label that can be updated below the dropdown menu
                        
                        self.wp_status = 'Pathway image may take several seconds to a minute to load...\n'
                        self.wp_status += '(images saved to "WikiPathways" folder in input directory)'
                        try: self.label_status_name.set(self.wp_status)
                        except Exception:
                            self.label_status_name = StringVar()
                            self.label_status_name.set(self.wp_status)
                            self.invokeStatusLabel() ### Invoke a new label indicating that the database is loading
                    except Exception:
                        None
                        
        if option == 'wp_id_selection':
            ### Reset any manually input WPID if a new pathway is selected from dropdown
            try: self.entry_field.setentry('') 
            except Exception: null=[]
            
    def viewPNGFile(self):
        """ View PNG file within a PMW Tkinter frame """
        import ImageTk ### HAVE TO CALL HERE TO TRIGGER AN ERROR - DON'T WANT THE TopLevel to open otherwise
        png_file_dir = self.graphic_link['WP']
        img = ImageTk.PhotoImage(file=png_file_dir)
        
        tl = Toplevel()
        sf = PmwFreeze.ScrolledFrame(tl, labelpos = 'n', label_text = '',
                usehullsize = 1, hull_width = 800, hull_height = 550)
        sf.pack(padx = 0, pady = 0, fill = 'both', expand = 1)
        frame = sf.interior()

        tl.title(png_file_dir)
        can = Canvas(frame)
        can.pack(fill=BOTH, padx = 0, pady = 0)
        w = img.width()
        h = height=img.height()
        
        can.config(width=w, height=h)        
        can.create_image(2, 2, image=img, anchor=NW)
        tl.mainloop()
        
    def openPNGImage(self):
        png_file_dir = self.graphic_link['WP']
        if os.name == 'nt':
            try: os.startfile('"'+png_file_dir+'"')
            except Exception:  os.system('open "'+png_file_dir+'"')
        elif 'darwin' in sys.platform: os.system('open "'+png_file_dir+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+png_file_dir+'"')   
        
    def __init__(self, parent, option_db, option_list): 
        self._parent = parent; self._option_list = option_list
        self._user_variables = user_variables
        self.default_dir = PathDir; self.default_file = PathFile
        self.option_db = option_db
        filename = 'Config/icon.gif'
        fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='top'); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        self.pathdb={}; use_scroll = 'no'

        label_text_str = "\nGO-Elite Program Options"
        if os.name == 'nt': height = 400; width = 420
        else: height = 350; width = 480
        if option_db == 'ViewWikiPathways': width = 520
        use_scroll = 'yes'
        if 'run_from_scratch' in option_list or 'modifyDBs1' in option_list:  
            if 'modifyDBs1' in option_list and os.name != 'nt': width = 460; height = 350
            else: width = 400; height = 325
        if 'selected_species3' in option_list:
            height -= 30
        if os.name != 'nt':height+=50; width+=50
        self.sf = PmwFreeze.ScrolledFrame(self._parent,
                labelpos = 'n', label_text = label_text_str,
                usehullsize = 1, hull_width = width, hull_height = height)
        self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
        self.frame = self.sf.interior()
        parent_type = self.frame
        self.parent_type = parent_type
        i = 0; object_directions = ['top','bottom','up','down']
        if option_db == 'ViewWikiPathways':
            self.ViewWikiPathways()
        for option in option_list:
          if option in self.option_db:
            od = self.option_db[option]; self.title = od.Display(); description = od.Description()      
            self.display_options = od.AnalysisOptions()
            if 'radio' in od.DisplayObject() and self.display_options != ['NA']:           
                ### Create and pack a RadioSelect widget, with radiobuttons.
                self._option = option
                def radiocallback(tag,callback=self.callback,option=option): callback(tag,option)
                radiobuttons = PmwFreeze.RadioSelect(parent_type,                       
                        buttontype = 'radiobutton', orient = 'vertical',
                        labelpos = 'w', command = radiocallback, label_text = self.title,
                        hull_borderwidth = 2, hull_relief = 'ridge',
                );
                if description in object_directions: direction = description ### can be used to store directions
                else: direction = 'top'
                radiobuttons.pack(side = direction, expand = 1, padx = 10, pady = 5)
                    
                ### print self.display_options
                ### Add some buttons to the radiobutton RadioSelect.
                for text in self.display_options:
                    if text != ['NA']: radiobuttons.add(text)
                self.default_option = od.DefaultOption()
                radiobuttons.invoke(self.default_option)
                
            if 'button' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                ### Create and pack a horizontal RadioSelect widget.
                self.default_option = od.DefaultOption()
                
                if mac_print_mode == 'yes': 
                    button_type = 'radiobutton'
                    if self._option == 'run_from_scratch':
                        self.title = self.title[:int(len(self.title)/1.8)] ### On a mac, the spaced text is too far out
                    if self._option == 'modifyDBs1':
                        self.title = self.title[:int(len(self.title)/1.3)] ### On a mac, the spaced text is too far out
                else: button_type = 'button'
                def buttoncallback(tag,callback=self.callback,option=option):
                    callback(tag,option)
                horiz = PmwFreeze.RadioSelect(parent_type, buttontype = button_type, orient = 'vertical',
                        labelpos = 'w', command = buttoncallback,
                        label_text = self.title, frame_borderwidth = 2,
                        frame_relief = 'ridge'
                ); horiz.pack(fill = 'x', padx = 10, pady = 10)

                ### Add some buttons to the horizontal RadioSelect
                for text in self.display_options:
                    if text != ['NA']: horiz.add(text)
                horiz.invoke(self.default_option)

            if ('folder' in od.DisplayObject() or 'file' in od.DisplayObject()) and self.display_options != ['NA']:
              proceed = 'yes'
              if option == 'mappfinder_dir' and run_mappfinder == 'yes': proceed = 'no'
              if proceed == 'yes':
                self._option = option

                if option == 'output_dir':
                    if run_mappfinder == 'no':
                        notes = ""#"note: if not selected, 'input/MAPPFinder' will\nbe used and results written to 'output/' "
                        self.title = string.replace(self.title,'output','pre-computed ORA')
                    else: notes = od.Description()
                else: notes = od.Description()
                
                group = PmwFreeze.Group(parent_type,tag_text = self.title)
                group.pack(fill = 'both', expand = 1, padx = 10, pady = 2)
                
                def filecallback(callback=self.callback,option=option): self.getPath(option)             
                entrytxt = StringVar(); #self.entrytxt.set(self.default_dir)
                default_option = string.replace(od.DefaultOption(),'---','')
                entrytxt.set(default_option)
                self.pathdb[option] = entrytxt
                self._user_variables[option] = default_option
                
                #l = Label(group.interior(), text=self.title); l.pack(side=LEFT)        
                entry = Entry(group.interior(),textvariable=self.pathdb[option]);
                entry.pack(side='left',fill = 'both', expand = 1, padx = 10, pady = 2)
                button = Button(group.interior(), text="select "+od.DisplayObject(), width = 10, fg="red", command=filecallback)
                button.pack(side=LEFT, padx = 2,pady = 2)                    

                #print option,run_mappfinder, self.title, self.default_option
                if len(notes)>0: ln = Label(parent_type, text=notes,fg="blue"); ln.pack(padx = 10)
                
            if 'drop-down' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                self.default_option = self.display_options
                def comp_callback1(tag,callback=self.callback,option=option):
                    callback(tag,option)
                    
                self.comp = PmwFreeze.OptionMenu(parent_type,
                    labelpos = 'w', label_text = self.title,                                                 
                    items = self.default_option, command = comp_callback1)
                
                if 'species' in option:
                    if 'selected_species2' in option:
                        self.speciescomp2 = self.comp; self.speciescomp2.pack(anchor = 'w', padx = 10, pady = 0)
                    elif 'selected_species3' in option:
                        self.speciescomp3 = self.comp; self.speciescomp3.pack(anchor = 'w', padx = 10, pady = 0)
                    else: self.speciescomp = self.comp; self.speciescomp.pack(anchor = 'w', padx = 10, pady = 0)
                    try: self.speciescomp.invoke(od.DefaultOption()) ###Just pick the first option
                    except Exception: self.speciescomp.invoke(self.default_option[0])
                else:
                    self.comp.pack(anchor = 'w', padx = 10, pady = 0)
                    try: self.comp.invoke(od.DefaultOption()) ###Just pick the first option
                    except Exception: self.comp.invoke(self.default_option[0])
                    if len(od.Description())>0 and od.Description() != 'top':
                        ln = Label(parent_type, text=od.Description(),fg="blue"); ln.pack(padx = 10)
                    if option == 'selected_version':
                        notes = 'Note:   Available species may vary based on database selection and\n'
                        notes+= '"Plus" versions have additional Affymetrix & EntrezGene relationships\n'
                        ln = Label(parent_type, text=notes,fg="blue"); ln.pack(padx = 10)

            if 'comboBox' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                self.default_option = self.display_options
                def comp_callback1(tag,callback=self.callbackComboBox,option=option):
                    callback(tag,option)

                self.comp = PmwFreeze.ComboBox(parent_type,
                    labelpos = 'w', dropdown=1, label_text = self.title,
                    unique = 0, history = 0,
                    scrolledlist_items = self.default_option,
                    selectioncommand = comp_callback1)
                        
                if 'species' in option:
                    if 'selected_species2' in option:
                        self.speciescomp2 = self.comp; self.speciescomp2.pack(anchor = 'w', padx = 10, pady = 0)
                        try: self.speciescomp2.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.speciescomp2.invoke())
                        except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                        try:
                            self.speciescomp2.selectitem(od.DefaultOption()) 
                            self.callbackComboBox(od.DefaultOption(),option) 
                        except Exception:
                            self.speciescomp2.selectitem(self.default_option[0]) ###Just pick the first option
                            self.callbackComboBox(self.default_option[0],option)
                    elif 'selected_species3' in option:
                        self.speciescomp3 = self.comp; self.speciescomp3.pack(anchor = 'w', padx = 10, pady = 0)
                        try: self.speciescomp3.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.speciescomp3.invoke())
                        except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                        try:
                            self.speciescomp3.selectitem(od.DefaultOption()) ###Just pick the first option
                            self.callbackComboBox(od.DefaultOption(),option) 
                        except Exception:
                            self.speciescomp3.selectitem(self.default_option[0])
                            self.callbackComboBox(self.default_option[0],option)
                    else:
                        self.speciescomp = self.comp; self.speciescomp.pack(anchor = 'w', padx = 10, pady = 0)
                        try: self.speciescomp.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.speciescomp.invoke())
                        except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                        try:
                            self.speciescomp.selectitem(od.DefaultOption())
                            self.callbackComboBox(od.DefaultOption(),option) 
                        except Exception:
                            self.speciescomp.selectitem(self.default_option[0])
                            self.callbackComboBox(self.default_option[0],option)
                else:
                    self.combo = self.comp ### has to be a unique combo box to refer to itself in the component call below
                    self.combo.pack(anchor = 'w', padx = 10, pady = 1)
                    try: self.combo.component('entryfield_entry').bind('<Button-1>', lambda event, self=self: self.combo.invoke())
                    except Exception: None ### Above is a slick way to force the entry field to be disabled and invoke the scrolledlist
                    try:
                        self.combo.selectitem(od.DefaultOption())
                        self.callbackComboBox(od.DefaultOption(),option) 
                    except Exception:
                        self.combo.selectitem(self.default_option[0])
                        self.callbackComboBox(self.default_option[0],option)
                    if len(od.Description())>0 and od.Description() != 'top':
                        ln = Label(parent_type, text=od.Description(),fg="blue"); ln.pack(padx = 10)
                    if option == 'selected_version':
                        notes = 'Note:   Available species may vary based on database selection and\n'
                        notes+= '"Plus" versions have additional Affymetrix & EntrezGene relationships\n'
                        ln = Label(parent_type, text=notes,fg="blue"); ln.pack(padx = 10)
                        
            if 'label' in od.DisplayObject() and self.display_options != ['NA']:
                ln = Label(parent_type, text=self.title,fg="blue"); ln.pack(padx = 10)
                        
            if 'enter' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                ### Create and pack a horizontal RadioSelect widget.
                self.default_option = od.DefaultOption()
                
                #print self.default_option, self.title; kill

                def custom_validate(tag,custom_validate=self.custom_validate,option=option):
                    validate = custom_validate(tag,option)
                def custom_validate_p(tag,custom_validate_p=self.custom_validate_p,option=option):
                    validate = custom_validate_p(tag,option)
                    #print [validate], tag, option
                try:
                    if float(self.default_option) <= 1: use_method = 'p'
                    else: use_method = 'i'
                except ValueError:
                   use_method = 'i'
                if use_method == 'p':
                    self.entry_field = PmwFreeze.EntryField(parent_type,
                            labelpos = 'w',
                            label_text = self.title,
                            validate = custom_validate, 
                            value = self.default_option, hull_borderwidth = 2, hull_relief = 'ridge'
                    ); self.entry_field.pack(fill = 'x', expand = 1, padx = 10, pady = 10)
                if use_method == 'i':
                    self.entry_field = PmwFreeze.EntryField(parent_type,
                            labelpos = 'w',
                            label_text = self.title,
                            validate = custom_validate,
                            value = self.default_option, hull_borderwidth = 2, hull_relief = 'ridge'
                    ); self.entry_field.pack(fill = 'x', expand = 1, padx = 10, pady = 10)
                    if len(od.Description())>0:
                        ln = Label(parent_type, text=od.Description(),fg="red"); ln.pack(padx = 10)
            if 'multiple-checkbox' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                self.default_option = od.DefaultOption()
                
                ### Create and pack a vertical RadioSelect widget, with checkbuttons.
                self.checkbuttons = PmwFreeze.RadioSelect(parent_type,
                        buttontype = 'checkbutton', orient = 'vertical',
                        labelpos = 'w', command = self.checkbuttoncallback,
                        label_text = self.title, hull_borderwidth = 2, hull_relief = 'ridge',
                ); self.checkbuttons.pack(side = 'top', expand = 1, padx = 10, pady = 10)

                ### Add some buttons to the checkbutton RadioSelect.
                for text in self.display_options:
                     if text != ['NA']: self.checkbuttons.add(text)
                self.checkbuttons.invoke(self.default_option)
                self.checkbuttons.invoke(self.default_option2)
                
            if 'single-checkbox' in od.DisplayObject() and self.display_options != ['NA']:
                self._option = option
                proceed = 'yes'
                """if option == 'export_splice_index_values':
                    if analysis_method != 'splicing-index': proceed = 'no' ### only export corrected constitutive ratios if splicing index method chosen"""
                if proceed == 'yes':
                    self.default_option = od.DefaultOption()
                    
                    if self.default_option != 'NA':
                        def checkbuttoncallback(tag,state,checkbuttoncallback=self.checkbuttoncallback,option=option):
                            checkbuttoncallback(tag,state,option)                
                        ### Create and pack a vertical RadioSelect widget, with checkbuttons. hull_relief
                        self.checkbuttons = PmwFreeze.RadioSelect(parent_type,
                                buttontype = 'checkbutton', command = checkbuttoncallback,
                                hull_borderwidth = 2,
                        ); self.checkbuttons.pack(side = 'top', expand = 1, padx = 10, pady = 10)

                        ### Add some buttons to the checkbutton RadioSelect.
                        self.checkbuttons.add(self.title)
                        if self.default_option == 'yes': self.checkbuttons.invoke(self.title)
                        else: self._user_variables[option] = 'no'
            i+=1 ####Keep track of index
        if len(option_list)>0:
            if 'process_go' in option_list: ### For the CEL file selection window, provide a link to get Annotation files
                button_text = 'Download Annotation CSVs'; d_url = 'http://www.affymetrix.com/support/technical/byproduct.affx?cat=arrays'
                self.d_url = d_url; text_button = Button(self._parent, text=button_text, command=self.Dlinkout); text_button.pack(side = 'left', padx = 5, pady = 5)
            if 'permutation' in option_list:
                #ScrolledCanvas().mainloop()
                systemCodes_win = Button(parent_type, text="GO-Elite Supported System Codes", command=self.systemCodes)
                systemCodes_win.pack(side = 'bottom', padx = 5, pady = 5)
            if 'new_species_name' in option_list:
                button_text = 'Lookup Species TaxID'; d_url = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=taxonomy'
                self.d_url = d_url; text_button = Button(self._parent, text=button_text, command=self.Dlinkout); text_button.pack(side = 'left', padx = 5, pady = 5)
                
            continue_to_next_win = Button(text = 'Continue', command = self._parent.destroy)
            continue_to_next_win.pack(side = 'right', padx = 10, pady = 10)
    
            back_button = Button(self._parent, text="Back", command=self.goBack) 
            back_button.pack(side = 'right', padx =10, pady = 5)
            
            quit_win = Button(self._parent, text="Quit", command=self.quit) 
            quit_win.pack(side = 'right', padx =10, pady = 5)
    
            button_text = 'Help'
            url = 'http://www.genmapp.org/go_elite/help_main.htm'; self.url = url
            pdf_help_file = 'Documentation/GO-Elite_Manual.pdf'; pdf_help_file = filepath(pdf_help_file); self.pdf_help_file = pdf_help_file
    
            try: help_button = Button(self._parent, text=button_text, command=self.GetHelpTopLevel); help_button.pack(side = 'left', padx = 5, pady = 5)
            except Exception: help_button = Button(self._parent, text=button_text, command=self.linkout); help_button.pack(side = 'left', padx = 5, pady = 5)
    
            self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
            self._parent.mainloop()

    def goBack(self):
        self._parent.destroy()
        if 'filter_method' in self._option_list: run_parameter = 'ORA',self._user_variables
        else: run_parameter = 'skip'
        reload(GO_Elite); GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
        
    def GetHelpTopLevel(self):
        message = ''
        self.message = message; self.online_help = 'Online Documentation'; self.pdf_help = 'Local PDF File'
        tl = Toplevel(); self._tl = tl; nulls = '\t\t\t\t'; tl.title('Please select one of the options')
        self.sf = PmwFreeze.ScrolledFrame(self._tl,
                labelpos = 'n', label_text = '',
                usehullsize = 1, hull_width = 220, hull_height = 150)
        self.sf.pack(padx = 10, pady = 10, fill = 'both', expand = 1)
        self.frame = self.sf.interior()
        group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Options')
        group.pack(fill = 'both', expand = 1, padx = 20, pady = 10)
        l1 = Label(group.interior(), text=nulls);  l1.pack(side = 'bottom')
        text_button2 = Button(group.interior(), text=self.online_help, command=self.openOnlineHelp); text_button2.pack(side = 'top', padx = 5, pady = 5) 
        try: text_button = Button(group.interior(), text=self.pdf_help, command=self.openPDFHelp); text_button.pack(side = 'top', padx = 5, pady = 5)
        except Exception: text_button = Button(group.interior(), text=self.pdf_help, command=self.openPDFHelp); text_button.pack(side = 'top', padx = 5, pady = 5)
        tl.mainloop()
    def openPDFHelp(self):
        if os.name == 'nt':
            try: os.startfile('"'+self.pdf_help_file+'"')
            except Exception:  os.system('open "'+self.pdf_help_file+'"')
        elif 'darwin' in sys.platform: os.system('open "'+self.pdf_help_file+'"')
        elif 'linux' in sys.platform: os.system('xdg-open "'+self.pdf_help_file+'"')   
        self._tl.destroy()
    def openOnlineHelp(self):
        try: webbrowser.open(self.url)
        except Exception: null=[]
        self._tl.destroy()
    
    def linkout(self):
        try: webbrowser.open(self.url)
        except Exception: null=[]
    def Dlinkout(self):
        try: webbrowser.open(self.d_url)
        except Exception: null=[]
            
    def setvscrollmode(self, tag):
        self.sf.configure(vscrollmode = tag)

    def systemCodes(self):
        selected_species = self._user_variables['species']
        species_systems = getSpeciesSystems(selected_species)
        GO_Elite.sourceData()
        about = "System Name"+'\t\t'+"System Code"+'\t\t'+"MOD"+'\n'
        system_list=[]
        for system_name in system_codes:
            sc = system_codes[system_name]
            if system_name in species_systems: ### Restrict by species
                system_list.append([system_name,sc.SystemCode(),sc.MOD()])
        system_list.sort()
        for (system_name,system_code,mod) in system_list:
            filler = '  '; filler2=''
            if os.name == 'nt':
                if platform.win32_ver()[0] != 'Vista': val = 9; val2 = 14; val3 = 19
                else: val = 12; val2 = 17; val3 = 20
            else: val = 12; val2 = 17; val3 = 20
            if os.name == 'nt':
                if len(system_name)<val: filler = val2-len(system_name); filler = filler*' '+'  '
            if len(mod)<val: filler2 = val3-len(mod); filler2 = filler2*' '
            about+= system_name+filler+'\t\t'+system_code+filler2+'\t\t'+mod+'\n'
        #tkMessageBox.showinfo("Available GO-Elite System Codes",about,parent=self._parent)

        dialog = PmwFreeze.TextDialog(self._parent, scrolledtext_labelpos = 'n',
                title = 'GO-Elite System Codes', defaultbutton = 0,
                label_text = 'Available System Information for '+selected_species)
        dialog.insert('end', about)
            
    def info(self):
        tkMessageBox.showinfo("title","message",parent=self._parent)

    def deleteWindow(self):
        #tkMessageBox.showwarning("Quit","Use 'Quit' button to end program!",parent=self._parent)
        self._parent.destroy(); sys.exit() ### just quit instead

    def quit(self):
        #print "quit starts"
        #print "cleaning up things..."
        self._parent.quit()
        self._parent.destroy()
        sys.exit()
        #print "quit ends"

    def chooseDirectory(self,option):
        tag = tkFileDialog.askdirectory(parent=self._parent)
        #print option,tag
        self._user_variables[option] = tag
        #tkFileDialog.askopenfile(parent=self._parent)

    def getPath(self,option):
        if 'dir' in option or 'folder' in option:
            try: dirPath = tkFileDialog.askdirectory(parent=self._parent,initialdir=self.default_dir)
            except Exception: 
                self.default_dir = ''
                try: dirPath = tkFileDialog.askdirectory(parent=self._parent,initialdir=self.default_dir)
                except Exception: 
                    try: dirPath = tkFileDialog.askdirectory(parent=self._parent)
                    except Exception: dirPath=''
            self.default_dir = dirPath
            entrytxt = self.pathdb[option]
            entrytxt.set(dirPath)
            self._user_variables[option] = dirPath
            try: file_location_defaults['PathDir'].SetLocation(dirPath)
            except Exception: null = None
            exportDefaultFileLocations(file_location_defaults)
            ### Allows the option_db to be updated for the next round (if an error is encountered)
            
        if 'file' in option:
            try: tag = tkFileDialog.askopenfile(parent=self._parent,initialdir=self.default_file)
            except Exception: 
                self.default_file = ''
                try: tag = tkFileDialog.askopenfile(parent=self._parent,initialdir=self.default_file)
                except Exception: 
                    try: tag = tkFileDialog.askopenfile(parent=self._parent)
                    except Exception: tag=''
            try: filePath = tag.name #initialdir=self.default_dir
            except AttributeError: filePath = ''
            filePath_dir = string.join(string.split(filePath,'/')[:-1],'/')
            self.default_file = filePath_dir
            entrytxt = self.pathdb[option]
            entrytxt.set(filePath)
            self._user_variables[option] = filePath
            try: file_location_defaults['PathFile'].SetLocation(filePath_dir)
            except Exception: null = None
            exportDefaultFileLocations(file_location_defaults)
    
    def Report(self,tag,option):
        output = tag
        return output
    def __repr__(self,tag,option): return self.Report(tag,option)
    
    def Results(self): return self._user_variables

    def custom_validate(self, text, option):
        #print [option],'text:', text
        self._user_variables[option] = text
        try:
            text = float(text);return 1
        except ValueError: return -1

    def custom_validate_p(self, text, option):
        #print [option],'text:', text
        self._user_variables[option] = text
        try:
            text = float(text)
            if text <1:return 1
            else:return -1
        except ValueError:return -1

    def callback(self, tag, option):
        #print 'Button',[option], tag,'was pressed.'
        self._user_variables[option] = tag
        if option == 'ORA_algorithm':
            if tag == 'Permute p-value':
                try: self.entry_field.setentry('2000')
                except Exception: null=[]
                self._user_variables['permutation'] = '2000'
            elif tag == 'Fisher Exact Test':
                try: self.entry_field.setentry('NA')
                except Exception: null=[]
                self._user_variables['permutation'] = '0'
        if option == 'dbase_version':
            ###Export new species info
            exportDBversion(tag)
            current_species_names = getSpeciesList()
            ### THIS TOOK FOR EVER TO FIND!!!!!! setitems of the PMW object resets the value list
            if 'permutation' in self._option_list:
                try: self.speciescomp.setitems(current_species_names)
                except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            else:
                try: self.speciescomp.setitems(['all-supported','New Species','---']+current_species_names)
                except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            for i in self._option_list:
                var = 'proceed'
                if 'species' in i: ### Necessary if the user changes dbase_version and selects continue to accept the displayed species name (since it's note directly invoked)
                    if 'species' in self._user_variables:
                        if self._user_variables[i] in current_species_names: var = None
                    if var == 'proceed':
                        try: self._user_variables[i] = current_species_names[0]
                        except Exception: null = []
        elif option == 'selected_version':
            current_species_names = db_versions[tag]
            current_species_names.sort()
            try: self.speciescomp.setitems(['---']+current_species_names)
            except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try: self.speciescomp2.setitems(['---']+current_species_names)
            except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try: self.speciescomp3.setitems(['---']+current_species_names)
            except Exception: null = [] ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated

    def callbackComboBox(self, tag, option):
        """ Similiar to the above, callback, but ComboBox uses unique methods """
        #print 'Button',[option], tag,'was pressed.'
        self._user_variables[option] = tag

        if option == 'selected_version':
            current_species_names = db_versions[tag]
            current_species_names.sort()
            current_species_names = ['---']+current_species_names
            species_option = current_species_names[0]
            try:
                self.speciescomp._list.setlist(current_species_names) ### This is the way we set a new list for ComboBox
                ### Select the best default option to display (keep existing or re-set)
                if 'selected_species1' in self._user_variables: ### If this is the species downloader
                    species_option = 'selected_species1'
                else:
                    for i in self._user_variables:
                        if 'species' in i: species_option = i
                default = self.getBestDefaultSelection(species_option,current_species_names)
                self.speciescomp.selectitem(default)
            except Exception: None ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try:
                self.speciescomp2._list.setlist(current_species_names)
                default = self.getBestDefaultSelection('selected_species2',current_species_names)
                self.speciescomp2.selectitem(default)
            except Exception: None ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
            try:
                self.speciescomp3._list.setlist(current_species_names)
                default = self.getBestDefaultSelection('selected_species3',current_species_names)
                self.speciescomp3.selectitem(default)
            except Exception: None ### Occurs before speciescomp is declared when dbase_version pulldown is first intiated
        
    def getBestDefaultSelection(self,option,option_list):
        default = option_list[0] ### set the default to the first option listed
        if option in self._user_variables:
            selected = self._user_variables[option]
            if selected in option_list: ### If selected species exists in the new selected version of EnsMart
                default = selected
            else:
                self._user_variables[option] = default ### Hence, the default has changed, so re-set it
        return default
        
    def multcallback(self, tag, state):
        if state: action = 'pressed.'
        else: action = 'released.'
        """print 'Button', tag, 'was', action, \
                'Selection:', self.multiple.getcurselection()"""
        self._user_variables[option] = tag

    def checkbuttoncallback(self, tag, state, option):
        if state: action = 'pressed.'
        else: action = 'released.'
        """print 'Button',[option], tag, 'was', action, \
                'Selection:', self.checkbuttons.getcurselection()"""
        if state==0: tag2 = 'no'
        else: tag2 = 'yes'
        #print '---blahh', [option], [tag], [state], [action], [self.checkbuttons.getcurselection()]
        self._user_variables[option] = tag2

def getSpeciesList():
    try: current_species_dirs = unique.read_directory('/Databases')
    except Exception: ### Occurs when the version file gets over-written with a bad directory name
        try:
            ### Remove the version file and wipe the species file
            os.remove(filepath('Config/version.txt'))
            raw = export.ExportFile('Config/species.txt'); raw.close()
            os.mkdir(filepath('Databases'))
        except Exception: null = []
        elite_db_versions = returnDirectoriesNoReplace('/Databases')
        try: exportDBversion(elite_db_versions[0])
        except Exception: exportDBversion('')
        current_species_dirs = unique.read_directory('/Databases')

    current_species_names=[]
    for species in species_codes:
        if species_codes[species].SpeciesCode() in current_species_dirs: current_species_names.append(species)
    current_species_names.sort()
    return current_species_names

def exportDBversion(db_version):
    today = str(datetime.date.today()); today = string.split(today,'-'); today = today[1]+'/'+today[2]+'/'+today[0]
    OBO_import.exportVersionData(db_version,today,'Config/')
    
def getSpeciesSystems(species):
    species_code = species_codes[species].SpeciesCode()
    import_dir1 = '/Databases/'+species_code+'/uid-gene'
    import_dir2 = '/Databases/'+species_code+'/gene-mapp'
    import_dir3 = '/Databases/'+species_code+'/gene-go'
    try: uid_gene_list = read_directory(import_dir1)
    except Exception: uid_gene_list=[]
    try: gene_mapp_list = read_directory(import_dir2)
    except Exception: gene_mapp_list=[]
    try: gene_go_list = read_directory(import_dir3)
    except Exception: gene_go_list=[]
    uid_gene_list+=gene_mapp_list+gene_go_list
    systems=[]
    for file in uid_gene_list:
        file = string.replace(file,'.txt','')
        systems += string.split(file,'-')
    systems = unique.unique(systems)
    return systems

class SpeciesData:
    def __init__(self, abrev, species, systems, taxid):
        self._abrev = abrev; self._species = species; self._systems = systems; self._taxid = taxid
    def SpeciesCode(self): return self._abrev
    def SpeciesName(self): return self._species
    def Systems(self): return self._systems
    def TaxID(self): return self._taxid
    def __repr__(self): return self.SpeciesCode()+'|'+self.SpeciesName()
        
def importSpeciesInfo():
    try:
        if integrate_online_species == 'yes': filename = 'Config/species_all.txt'
        else: filename = 'Config/species.txt'
    except Exception: filename = 'Config/species.txt'
    fn=filepath(filename); global species_list; species_list=[]; global species_codes; species_codes={}; x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        try:
            abrev,species,taxid,compatible_mods = string.split(data,'\t')
        except Exception:
            if '!DOCTYPE': print_out = "A internet connection could not be established.\nPlease fix the problem before proceeding."
            else: print_out = "Unknown file error encountered."
            raw = export.ExportFile(fn); raw.close(); GO_Elite.importGOEliteParameters('skip'); sys.exit()
        if x==0: x=1
        else:
            compatible_mods = string.split(compatible_mods,'|')
            species_list.append(species)
            sd = SpeciesData(abrev,species,compatible_mods,taxid)
            species_codes[species] = sd
    return species_codes

def remoteSpeciesInfo(integrate_online):
    global integrate_online_species; integrate_online_species = integrate_online; species_names={}
    importSpeciesInfo()
    if len(species_codes) == 0:
        integrate_online_species = 'yes'
        importSpeciesInfo() ### Occurs when an unknown error erases the species file
        exportSpeciesInfo(species_codes) ### Re-sets this database
    for species in species_codes:
        sd = species_codes[species]
        species_names[sd.SpeciesCode()] = sd
    return species_names
    
def exportSpeciesInfo(species_codes):
    fn=filepath('Config/species.txt'); data = open(fn,'w'); x=0
    header = string.join(['species_code','species_name','tax_id','compatible_algorithms'],'\t')+'\n'
    data.write(header)
    for species in species_codes:
        if 'New Species' not in species and 'all-' not in species and species != '':
            sd = species_codes[species]
            mods = string.join(sd.Systems(),'|')
            values = [sd.SpeciesCode(),sd.SpeciesName(),sd.TaxID(),mods]
            values = string.join(values,'\t')+'\n'
            data.write(values)
    data.close()

def exportArrayVersionInfo(species_array_db):
    fn=filepath('Config/array_versions.txt'); data = open(fn,'w'); x=0
    print 'Exporting:',fn
    for (species,db_version) in species_array_db:
        supported_arrays = unique.unique(species_array_db[(species,db_version)])
        supported_arrays = string.join(supported_arrays,'|')
        values = [species,db_version,supported_arrays]
        values = string.join(values,'\t')+'\n'
        data.write(values)
    data.close()

def exportSpeciesVersionInfo(species_archive_db):
    fn=filepath('Config/versions.txt'); data = open(fn,'w'); x=0
    print 'Exporting:',fn
    for species in species_archive_db:
        db_versions = species_archive_db[species]
        db_versions = string.join(db_versions,'|')
        values = [species,db_versions]
        values = string.join(values,'\t')+'\n'
        data.write(values)
    data.close()
    
def importOnlineDatabaseVersions():
    filename = 'Config/versions.txt'
    fn=filepath(filename); global db_versions; db_versions={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        try:
            species,versions = string.split(data,'\t')
            versions = string.split(versions,'|')
            for version in versions:
                try: db_versions[version].append(species)
                except KeyError: db_versions[version] = [species]
        except Exception: print [data]
    return db_versions

class SystemData:
    def __init__(self, syscode, sysname, mod):
        self._syscode = syscode; self._sysname = sysname; self._mod = mod
    def SystemCode(self): return self._syscode
    def SystemName(self): return self._sysname
    def MOD(self): return self._mod
    def __repr__(self): return self.SystemCode()+'|'+self.SystemName()+'|'+self.MOD()

def importSystemInfo():
    filename = 'Config/source_data.txt'; x=0
    fn=filepath(filename); mod_list=[]
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        if '!DOCTYPE' in data:
            fn2 = string.replace(fn,'.txt','_archive.txt')
            import shutil; shutil.copyfile(fn2,fn) ### Bad file was downloaded (with warning)
            importSystemInfo(); break
        else:
            try: sysname=t[0];syscode=t[1]
            except Exception: sysname=''
        try: mod = t[2]
        except Exception: mod = ''
        if x==0: x=1
        elif sysname != '':
            system_list.append(sysname)
            ad = SystemData(syscode,sysname,mod)
            if len(mod)>1: mod_list.append(sysname)
            system_codes[sysname] = ad
    return system_codes,system_list,mod_list

def remoteSystemInfo():
    global system_list; system_list=[]; global system_codes; system_codes={}
    system_codes,system_list,mod_list = importSystemInfo()
    return system_codes,system_list,mod_list

def exportSystemInfo():
    if len(system_codes)>0:
        filename = 'Config/source_data.txt'
        fn=filepath(filename); data = open(fn,'w')
        header = string.join(['System','SystemCode','MOD_status'],'\t')+'\n'
        data.write(header)
        for sysname in system_codes:
            ad = system_codes[sysname]
            values = string.join([sysname,ad.SystemCode(),ad.MOD()],'\t')+'\n'
            data.write(values)
        data.close()

def exportSystemInfoRemote(system_code_db):
    global system_codes; system_codes = system_code_db
    exportSystemInfo()
    
class FileLocationData:
    def __init__(self, status, location, species):
        self._status = status; self._location = location; self._species = species
    def Status(self): return self._status
    def Location(self): return self._location
    def SetLocation(self,location): self._location = location
    def Species(self): return self._species
    def __repr__(self): return self.Report()
    
def importDefaultFileLocations():
    filename = 'Config/default-files.csv'
    fn=filepath(filename); file_location_defaults={}
    for line in open(fn,'rU').readlines():
        line = string.replace(line,',','\t') ### Make tab-delimited (had to make CSV since Excel would impoperly parse otherwise)
        data = cleanUpLine(line)
        app,status,location,species = string.split(data,'\t')
        fl = FileLocationData(status, location, species)
        if species == 'all': file_location_defaults[app] = fl
        else:
            try: file_location_defaults[app].append(fl)
            except KeyError: file_location_defaults[app] = [fl]
    return file_location_defaults

def exportDefaultFileLocations(file_location_defaults):
    ### If the user supplies new defaults, over-write the existing
    fn=filepath('Config/default-files.csv'); data = open(fn,'w')
    for app in file_location_defaults:
        fl_list = file_location_defaults[app]
        try:
            for fl in fl_list:
                values = [app,fl.Status(),fl.Location(),fl.Species()]
                values = '"'+string.join(values,'","')+'"'+'\n'
                data.write(values)
        except Exception:
            fl = fl_list
            values = [app,fl.Status(),fl.Location(),fl.Species()]
            values = '"'+string.join(values,'","')+'"'+'\n'
            data.write(values)
    data.close()
    
class OptionData:
    def __init__(self,option,displayed_title,display_object,description,analysis_options,defaults):
        self._option = option; self._displayed_title = displayed_title; self._description = description
        self._analysis_options = analysis_options; self._display_object = display_object; self._defaults = defaults
    def Option(self): return self._option
    def Display(self): return self._displayed_title
    def DisplayObject(self): return self._display_object
    def Description(self): return self._description
    def AnalysisOptions(self): return self._analysis_options
    def DefaultOption(self): return self._defaults
    def setAnalysisOptions(self,analysis_options): self._analysis_options = analysis_options
    def setDefaultOption(self,defaults): self._defaults = defaults
    def __repr__(self): return self.Report()

def importUserOptions(type):
    filename = 'Config/options.txt'; option_db={}; option_group_db={}
    fn=filepath(filename); x=0
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        data = string.replace(data,'\k','\n') ###Used \k in the file instead of \n, since these are removed above
        t = string.split(data,'\t')
        option,displayed_title,display_object,group,description,analysis_options,defaults = t
        if x == 0:
           x = 1
        else:
            analysis_options = string.split(analysis_options,'|')
            od = OptionData(option,displayed_title,display_object,description,analysis_options,defaults)
            option_db[option] = od
            #option_list.append(option)
            try: option_group_db[group].append(option)
            except KeyError: option_group_db[group] = [option]
    return option_group_db,option_db 

class IndicatorWindow:
    def __init__(self,message,button_text):
        self.message = message; self.button_text = button_text
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')

        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()  
        quit_button = Button(parent, text='Quit', command=self.quit); quit_button.pack(side = 'bottom', padx = 5, pady = 5)
        text_button = Button(parent, text=self.button_text, command=parent.destroy); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        parent.mainloop()
        
    def quit(self): self._parent.quit(); self._parent.destroy(); sys.exit()
    def deleteWindow(self):
        tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)

class FeedbackWindow:
    def __init__(self,message,button_text,button_text2):
        self.message = message; self.button_text = button_text; self.button_text2 = button_text2
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Attention!!!')
        self._user_variables={}
        
        filename = 'Config/warning_big.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()
        text_button = Button(parent, text=self.button_text, command=self.button1); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        if button_text2 != '':
            text_button2 = Button(parent, text=self.button_text2, command=self.button2); text_button2.pack(side = 'bottom', padx = 5, pady = 5)
        parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        parent.mainloop()
        
    def button1(self): self._user_variables['button']=self.button_text; self._parent.destroy()
    def button2(self): self._user_variables['button']=self.button_text2; self._parent.destroy()
    def ButtonSelection(self): return self._user_variables
    def deleteWindow(self):
        tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
    
class ProcessCompletedWindow:
    def __init__(self,message,button_text):
        self.message = message; self.button_text = button_text
        parent = Tk(); self._parent = parent; nulls = '\t\t\t\t\t\t\t'; parent.title('Process Completed')

        filename = 'Config/icon.gif'; fn=filepath(filename); img = PhotoImage(file=fn)
        can = Canvas(parent); can.pack(side='left',padx = 10); can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)
        
        Label(parent, text='\n'+self.message+'\n'+nulls).pack()  
        quit_button = Button(parent, text='Quit', command=self.quit); quit_button.pack(side = 'bottom', padx = 5, pady = 5)
        text_button = Button(parent, text=self.button_text, command=parent.destroy); text_button.pack(side = 'bottom', padx = 5, pady = 5)
        parent.mainloop()
    def quit(self): self._parent.quit(); self._parent.destroy(); sys.exit()
    
class WarningWindow:
    def __init__(self,warning,window_name):
        tkMessageBox.showwarning(window_name, warning)

class InfoWindow:
    def __init__(self,dialogue,header):
        tkMessageBox.showinfo(header, dialogue)

class StringVarFile:
    def __init__(self,stringVar,window):
        self.__newline = 0; self.__stringvar = stringVar; self.__window = window
    def write(self,s):
        new = self.__stringvar.get()
        for c in s:
            #if c == '\n': self.__newline = 1
            if c == '\k': self.__newline = 1### This should not be found and thus results in a continous feed rather than replacing a single line
            else:
                if self.__newline: new = ""; self.__newline = 0
                new = new+c
        try: self.set(new)
        except Exception: None ### Not sure why this occurs
    def set(self,s): self.__stringvar.set(s); self.__window.update()   
    def get(self): return self.__stringvar.get()
                
class StatusWindow:
    def __init__(self,root,analysis,values):
        try:
            if debug_mode == 'yes':
                  if analysis == 'getAdditionalOnlineResources':
                      species_code,additional_resources = values
                      getAdditionalOnlineResources(species_code,additional_resources,root)
                  if analysis == 'EntrezGOExport':
                      tax_id,species_code,status,option_db,option_list,overwrite_entrezgo,rewrite_existing_EG = values
                      exportEntrezGO(tax_id,species_code,status,option_db,option_list,overwrite_entrezgo,rewrite_existing_EG,root)
                  if analysis == 'AffyCSV':
                      system_codes,species_code,species_full,incorporate_previous_associations,process_go,parse_genesets,integrate_affy_associations,overwrite_affycsv,get_mapped_file_version = values
                      extractAffymetrixCSVAnnotations(system_codes,species_code,species_full,incorporate_previous_associations,process_go,parse_genesets,integrate_affy_associations,overwrite_affycsv,get_mapped_file_version.root)
                  if analysis == 'UpdateOBO':
                      file_location_defaults,update_OBO,OBO_url = values
                      updateOBOfiles(file_location_defaults,update_OBO,OBO_url,root)
                  if analysis == 'getOnlineEliteConfig':
                      file_location_defaults = values
                      getOnlineEliteConfig(file_location_defaults,root)
                  if analysis == 'getOnlineEliteDatabase':
                      file_location_defaults,db_version,new_species_codes,download_obo,additional_resources = values
                      getOnlineEliteDatabase(file_location_defaults,db_version,new_species_codes,download_obo,additional_resources,root)
                  if analysis == 'DownloadEntrezGO':
                      file,dir,file_type,tax_id,species_code,overwrite_entrezgo,rewrite_existing_EG = values 
                      downloadFiles(file,dir,file_type,tax_id,species_code,overwrite_entrezgo,rewrite_existing_EG,root)
                  if analysis == 'EnsemblSQLImport':
                      species_code,species,child_dirs,externalDBName_list,overwrite_ensembl,rewrite_existing,force,ensembl_version,external_system = values
                      importEnsemblSQL(species_code,species,child_dirs,externalDBName_list,overwrite_ensembl,rewrite_existing,force,ensembl_version,external_system,root)
                  if analysis == 'copy':
                      file1,file2 = values
                      copyFiles(file1,file2,root)
                  if analysis == 'updateRelationshipFiles':
                      update_relationship_file,relationship_to_update,overwrite_relationships,species_code = values
                      addNewCustomRelationships(update_relationship_file,relationship_to_update,overwrite_relationships,species_code,root)
                  if analysis == 'updateAnnotationFiles':
                      update_annotation_file,annotations_to_update,overwrite_annotations,species_code = values
                      addNewCustomAnnotations(update_annotation_file,annotations_to_update,overwrite_annotations,species_code,root)
                  if analysis == 'getCurrentEnsemblSpecies':
                      getCurrentEnsemblSpecies(root)
                  if analysis == 'getVersionedEnsExternalDB':
                      species_full,ensembl_version = values
                      getVersionedEnsExternalDB(species_full,ensembl_version,root)
            else:
                  self._parent = root
                  root.title('GO-Elite 1.2.5 - Status Window')
                  statusVar = StringVar() ### Class method for Tkinter. Description: "Value holder for strings variables."
                  if os.name == 'nt': width = 575; height = 450
                  else: width = 650; height = 500
                  self.sf = PmwFreeze.ScrolledFrame(self._parent,
                        labelpos = 'n', label_text = 'Results Status Window',
                        usehullsize = 1, hull_width = width, hull_height = height)
                  self.sf.pack(padx = 5, pady = 1, fill = 'both', expand = 1)
                  self.frame = self.sf.interior()

                  group = PmwFreeze.Group(self.sf.interior(),tag_text = 'Output')
                  group.pack(fill = 'both', expand = 1, padx = 10, pady = 0)
                    
                  Label(group.interior(),width=180,height=452,justify=LEFT, bg='black', fg = 'white',anchor=NW,padx = 5,pady = 5, textvariable=statusVar).pack(fill=X,expand=Y)
                    
                  status = StringVarFile(statusVar,root) ###  Captures the stdout (or print) to the GUI instead of to the terminal
                  original_sys_out = sys.stdout ### Save the original stdout mechanism
                  sys.stdout = status
                  if analysis == 'getAdditionalOnlineResources':
                      species_code,additional_resources = values
                      root.after(100, getAdditionalOnlineResources(species_code,additional_resources,root))
                  if analysis == 'EntrezGOExport':
                      tax_id,species_code,status,option_db,option_list,overwrite_entrezgo,rewrite_existing_EG = values
                      root.after(100, exportEntrezGO(tax_id,species_code,status,option_db,option_list,overwrite_entrezgo,rewrite_existing_EG,root))
                  if analysis == 'AffyCSV':
                      system_codes,species_code,species_full,incorporate_previous_associations,process_go,parse_genesets,integrate_affy_associations,overwrite_affycsv,get_mapped_file_version = values
                      root.after(100, extractAffymetrixCSVAnnotations(system_codes,species_code,species_full,incorporate_previous_associations,process_go,parse_genesets,integrate_affy_associations,overwrite_affycsv,get_mapped_file_version,root))
                  if analysis == 'UpdateOBO':
                      file_location_defaults,update_OBO,OBO_url = values
                      root.after(100, updateOBOfiles(file_location_defaults,update_OBO,OBO_url,root))
                  if analysis == 'getOnlineEliteConfig':
                      file_location_defaults = values
                      root.after(100, getOnlineEliteConfig(file_location_defaults,root))
                  if analysis == 'getOnlineEliteDatabase':
                      file_location_defaults,db_version,new_species_codes,download_obo,additional_resources = values
                      root.after(100, getOnlineEliteDatabase(file_location_defaults,db_version,new_species_codes,download_obo,additional_resources,root))
                  if analysis == 'DownloadEntrezGO':
                      file,dir,file_type,tax_id,species_code,overwrite_entrezgo,rewrite_existing_EG = values 
                      root.after(100, downloadFiles(file,dir,file_type,tax_id,species_code,overwrite_entrezgo,rewrite_existing_EG,root))
                  if analysis == 'EnsemblSQLImport':
                      species_code,species,child_dirs,externalDBName_list,overwrite_ensembl,rewrite_existing,force,ensembl_version,external_system = values
                      importEnsemblSQL(species_code,species,child_dirs,externalDBName_list,overwrite_ensembl,rewrite_existing,force,ensembl_version,external_system,root)
                  if analysis == 'copy':
                      file1,file2 = values
                      root.after(100,copyFiles(file1,file2,root))
                  if analysis == 'updateRelationshipFiles':
                      update_relationship_file,relationship_to_update,overwrite_relationships,species_code = values
                      root.after(100,addNewCustomRelationships(update_relationship_file,relationship_to_update,overwrite_relationships,species_code,root))
                  if analysis == 'updateAnnotationFiles':
                      update_annotation_file,annotations_to_update,overwrite_annotations,species_code = values
                      root.after(100,addNewCustomAnnotations(update_annotation_file,annotations_to_update,overwrite_annotations,species_code,root))
                  if analysis == 'getCurrentEnsemblSpecies':
                      root.after(100,getCurrentEnsemblSpecies(root))
                  if analysis == 'getVersionedEnsExternalDB':
                      species_full,ensembl_version = values
                      root.after(100,getVersionedEnsExternalDB(species_full,ensembl_version,root))
                  #self._parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
                  self._parent.mainloop()
                  sys.stdout = original_sys_out ### Set this back to not capture print statements to the GUI anymore (if you don't no printing to the terminal)
                  try: self._parent.destroy()
                  except Exception: null = []
        except Exception,e:
            try:
                print traceback.format_exc()
                print_out = "Unknown error encountered during data processing.\nIf this error occurs again, please report to genmapp@gladstone.ucsf.edu."
                try: WarningWindow(print_out,'Error Encountered!'); self._parent.destroy()
                except Exception: print print_out
            except Exception: sys.exit()
        
    def deleteWindow(self): tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
    def quit(self): self._parent.quit(); self._parent.destroy(); sys.exit()
   
def addNewCustomRelationships(update_relationship_file,relationship_to_update,overwrite_relationships,species_code,root):
    try: status = gene_associations.addNewCustomRelationships(update_relationship_file,relationship_to_update,overwrite_relationships,species_code)
    except Exception, e: status = e
    if status == 'exported':
        print_out = "Imported relationships added to: "+relationship_to_update
        InfoWindow(print_out,'Process Completed')
    else:
        print status
        print_out = "Unknown file error encountered. Ensure the input file only has two columns."
        WarningWindow(print_out,'Error Encountered!')
    root.destroy()
    GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()

def addNewCustomAnnotations(update_annotation_file,annotations_to_update,overwrite_annotations,species_code,root):
    try: status = gene_associations.addNewCustomSystem(update_annotation_file,annotations_to_update,overwrite_annotations,species_code)
    except Exception, e: status = e; print e
    if status == 'exported':
        print_out = "Imported annotations added to "+species_code+' '+annotations_to_update
        InfoWindow(print_out,'Process Completed')
    else:
        print status
        print_out = "Unknown file error encountered. Ensure the input file only has three columns."
        WarningWindow(print_out,'Error Encountered!')
    root.destroy()
    GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()

def getCurrentEnsemblSpecies(root):
    global child_dirs; global ensembl_species; global ensembl_versions
    import EnsemblSQL
    try:
        print "Getting a list of current Ensembl species...one moment please"
        version = 'current'
        child_dirs, ensembl_species, ensembl_versions = EnsemblSQL.getCurrentEnsemblSpecies(version)
        ensembl_versions.reverse()
        root.destroy()
    except Exception:
        print_out = "A internet connection could not be established.\nPlease fix the problem before proceeding."
        WarningWindow(print_out,'Continue'); root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()

def filterExternalDBs(all_external_ids,externalDBName_list,external_ids,array_db):
    filtered_external_list=[]
    for name in externalDBName_list:
        if name in external_ids:
            id = external_ids[name]
            if id in all_external_ids:
                if name != 'GO': filtered_external_list.append(name)
    for array in array_db:
        if '\\N_' not in array: filtered_external_list.append(array)
    return filtered_external_list

def getVersionedEnsExternalDB(species_full,ensembl_version,root):
    try:
        import EnsemblSQL
        print "Getting a list of current Ensembl external databases...one moment please"
        child_dirs, ensembl_species, ensembl_versions = EnsemblSQL.getCurrentEnsemblSpecies(ensembl_version)
        try:
            ensembl_sql_dir,ensembl_sql_description_dir = child_dirs[species_full]
            ### Download the latest version of Ensembl
            try:
                EnsemblSQL.updateFiles(ensembl_sql_dir,'Config/','external_db.txt','yes')
                try: EnsemblSQL.updateFiles(string.replace(ensembl_sql_dir,'core','funcgen'),'Config/','array.txt','yes')
                except Exception:
                    raw = export.ExportFile('Config/array.txt'); raw.close()
                    print "No array relationships avaiable for",species_full
                root.destroy()
            except Exception:
                print_out = "Ensembl external database file not found."
                WarningWindow(print_out,'Continue'); root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
        except Exception:
            ### species not supported by Ensembl version
            print_out = "This species is not available for this version of\nEnsembl. Please try another version.."
            WarningWindow(print_out,'Continue'); root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
        
    except Exception:
        print_out = "A internet connection could not be established.\nPlease fix the problem before proceeding."
        WarningWindow(print_out,'Continue'); root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()

def copyFiles(file1,file2,root):
    print 'Copying file from:\n',file1
    print 'To:\n',file2
    data = export.ExportFile(file2) ### Ensures the directory exists
    shutil.copyfile(file1,file2)
    root.destroy()

def TimeStamp():
    time_stamp = time.localtime()
    year = str(time_stamp[0]); month = str(time_stamp[1]); day = str(time_stamp[2])
    if len(month)<2: month = '0'+month
    if len(day)<2: day = '0'+day
    return year+month+day

def deleteWPFiles():
    status = export.deleteFolder('BuildDBs/wikipathways')
    os.mkdir(filepath('BuildDBs/wikipathways'))
    
def extractAffymetrixCSVAnnotations(system_codes,species_codees,species_fulls,incorporate_previous_associations,process_go,parse_genesets,integrate_affy_associations,overwrite_affycsv,get_mapped_file_version,root):
    run_parameter = "Create/Modify Databases"
    
    if 'over-write previous' in overwrite_affycsv: overwrite_affycsv = 'over-write previous' ### The config name is longer
    if parse_genesets == 'no':
        for species_code in species_codees:
            print incorporate_previous_associations
            BuildAffymetrixAssociations.buildAffymetrixCSVAnnotations(species_code,incorporate_previous_associations,process_go,'no',integrate_affy_associations,overwrite_affycsv)
        print_out = 'Finished parsing the latest Affymetrix CSV annotations.'
        InfoWindow(print_out,'Update Complete!')
        continue_to_next_win = Button(text = 'Continue', command = root.destroy)
        continue_to_next_win.pack(side = 'right', padx = 10, pady = 10); root.mainloop()
        GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
    else:
        ### Coded a little weird, but allows the user to select mapped (Ensembl-inferred) versus non-mapped
        relationship_types = ['native','mapped']; ri=0
        if get_mapped_file_version == 'yes':
             relationship_types = ['mapped','mapped']

        ### Used when building a flat file from GPML zip file
        import gene_associations; all_species = 'no'
        try:
            gene_associations.convertAllGPML(species_codees,species_fulls) ### Downloads GPMLs and builds flat files
        except Exception:
            status = 'Unable to connect to http://www.wikipathways.org'
            try: WarningWindow(status,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
            except Exception: print status; sys.exit()
                        
        for relationship_type in relationship_types:
            index=0; ri+=1
            for species_code in species_codees:
                species_full = species_fulls[index]
                counts = BuildAffymetrixAssociations.importWikipathways(system_codes,incorporate_previous_associations,process_go,species_full,species_code,integrate_affy_associations,relationship_type,overwrite_affycsv)
                index+=1

        if counts == 0: print_out = 'No Affymetrix annotation files found, thus results based on existing meta file relationships.'; WarningWindow(print_out,'Update Incomplete!')
        else: print_out = 'Finished parsing the latest Wikipathways and Affymetrix annotations.'; InfoWindow(print_out,'Update Complete!')
        continue_to_next_win = Button(text = 'Continue', command = root.destroy)
        continue_to_next_win.pack(side = 'right', padx = 10, pady = 10); root.mainloop()
        GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
        
    print_out = 'Finished extracting Affymetrix probeset-gene associations file.'
    InfoWindow(print_out,'Update Complete!')  
    print 'Update Complete!'
    continue_to_next_win = Button(text = 'Continue', command = root.destroy)
    continue_to_next_win.pack(side = 'right', padx = 10, pady = 10); root.mainloop()
    GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
    
class MainMenu:
    def __init__(self):
        parent = Tk()
        self._parent = parent
        parent.title('GO-Elite: Introduction')
        self._user_variables={}
        filename = 'Config/logo.gif'
        fn=filepath(filename)
        img = PhotoImage(file=fn)
        can = Canvas(parent)
        can.pack(fill=BOTH)
        can.config(width=img.width(), height=img.height())        
        can.create_image(2, 2, image=img, anchor=NW)

        """      
        ### Create and pack a horizontal RadioSelect widget.
        def buttoncallback(tag,callback=self.callback):
            callback(tag)
        horiz = PmwFreeze.RadioSelect(parent,
                labelpos = 'w', command = buttoncallback,
                label_text = 'GO-Elite version 1.2.5 Main', frame_borderwidth = 2,
                frame_relief = 'ridge'
        ); horiz.pack(fill = 'x', padx = 10, pady = 10)
        for text in ['Continue']: horiz.add(text)
        """
        ### Add some buttons to the horizontal RadioSelect
        continue_to_next_win = Tkinter.Button(text = 'Begin Analysis', command = parent.destroy)
        continue_to_next_win.pack(side = 'bottom', padx = 5, pady = 5);

        info_win = Button(self._parent, text="About GO-Elite", command=self.info)
        info_win.pack(side = 'bottom', padx = 5, pady = 5)

        parent.protocol("WM_DELETE_WINDOW", self.deleteWindow)
        parent.mainloop()

    def info(self):
        about = 'GO-Elite version 1.2.5.\n'
        about+= 'GO-Elite is an open-source, freely available application covered under the\n'
        about+= 'Apache open-source license. Additional information can be found at:\n'
        about+= 'http://www.genmapp.org/go_elite\n'
        about+= '\nDeveloped by:\n\tNathan Salomonis\n\tBruce Conklin\nGladstone Institutes 2005-2012'
        tkMessageBox.showinfo("About GO-Elite",about,parent=self._parent)

    def deleteWindow(self):
        #tkMessageBox.showwarning("Quit Selected","Use 'Quit' button to end program!",parent=self._parent)
        self._parent.destroy(); sys.exit()

    def callback(self, tag):
        #print 'Button',[option], tag,'was pressed.'
        self._user_variables['continue'] = tag     
            
def addOnlineSpeciesDatabases():
    root = Tk(); analysis = 'getOnlineEliteConfig'
    values = file_location_defaults
    StatusWindow(root,analysis,values)

    importSystemInfo(); exportSystemInfo() ### By re-importing we incorporate new source data from the downloaded file    
    existing_species_codes = species_codes
    importSpeciesInfo()
    try: resource_list = importResourceList()
    except Exception: resource_list=['Failed']
    if len(species_codes) == 0:
        integrate_online_species = 'yes'
        importSpeciesInfo() ### Occurs when an unknown error erases the species file
    online_species = ['']
    for species in species_codes: online_species.append(species)
    online_species.sort()
    
    importOnlineDatabaseVersions(); db_version_list=[]
    for version in db_versions: db_version_list.append(version)
    db_version_list.sort(); db_version_list.reverse(); select_version = db_version_list[0]
    db_versions[select_version].sort()
    option_db['selected_species1'].setAnalysisOptions(['---']+db_versions[select_version])
    option_db['selected_species2'].setAnalysisOptions(['---']+db_versions[select_version])
    option_db['selected_species3'].setAnalysisOptions(['---']+db_versions[select_version])
    option_db['selected_version'].setAnalysisOptions(db_version_list)
    option_db['additional_resources'].setAnalysisOptions(['---']+resource_list)
    proceed = 'no'
    
    while proceed == 'no':
        root = Tk(); root.title('GO-Elite: Species Databases Available for Download')
        gu = GUI(root,option_db,option_group_db['OnlineDatabases'])
        db_version = gu.Results()['selected_version']
        exportDBversion(db_version)
        try: species1 = gu.Results()['selected_species1']
        except Exception: species1='---'
        try: species2 = gu.Results()['selected_species2']
        except Exception: species2='---'
        try: species3 = gu.Results()['selected_species3']
        except Exception: species3='---'
        try: download_obo = gu.Results()['download_obo']
        except Exception: download_obo='---'
        additional_resources = gu.Results()['additional_resources']
        new_species_list = [species1,species2,species3]; new_species_codes={}
        for species in new_species_list:
            if '---' not in species:
                try:
                    sc = species_codes[species].SpeciesCode();
                    new_species_codes[sc]=[]
                    existing_species_codes[species] = species_codes[species]
                except Exception: null=[]
        if len(new_species_codes) > 0 or download_obo == 'yes':
            root = Tk(); analysis = 'getOnlineEliteDatabase'
            values = file_location_defaults,db_version,new_species_codes,download_obo,additional_resources
            StatusWindow(root,analysis,values)
            proceed = 'yes'
        else:
            print_out = "Please select a species before continuing."
            IndicatorWindow(print_out,'Try Again')
    exportSpeciesInfo(existing_species_codes)
    integrate_online_species = 'no'
    
def getUserParameters(run_parameter):
    global root; global GO_Elite
    import GO_Elite
    if run_parameter == 'intro':
        try: MainMenu()
        except Exception:
            print_out = "\nCritical error encountered!!! This machine does not have either:\n"
            print_out += "1) Have the required Tcl/Tk components installed.\n"
            print_out += "2) Is being run from a compiled version that has critical incompatibilities your OS or hardware or\n"
            print_out += "3) Is being run from source-code in the same-directory as executable code resulting in a conflict\n"
            print_out += "\nIf any of these apply, we recommend downloading the Python source-code version of GO-Elite "
            print_out += "(installing necessary dependencies - see our Wiki or Documentation) or should be run from the web version of GO-Elite."
            print_out += "Otherwise, please contact GO-Elite support (genmapp@gladstone.ucsf.edu).\n\n"
            print_out += "Installation Wiki: http://code.google.com/p/go-elite/wiki/Installation\n\n"
            print print_out
        
            try:
                ### Create a log report of this
                try: log_file = filepath('GO-Elite_error-report.log')
                except Exception: log_file = filepath('/GO-Elite_error-report.log')
                log_report = open(log_file,'w');
                log_report.write(print_out)
                log_report.write(traceback.format_exc())
                log_report.close()
                ### Open this file
                if os.name == 'nt':
                    try: os.startfile('"'+log_file+'"')
                    except Exception:  os.system('open "'+log_file+'"')
                elif 'darwin' in sys.platform: os.system('open "'+log_file+'"')
                elif 'linux' in sys.platform: os.system('xdg-open "'+log_file+'/"')   
            except Exception: None
            sys.exit()
    
    
    global species; species=''; global user_variables; user_variables={}; global analysis_method; global array_type 
    ### Get default options for GO-Elite
    global run_mappfinder; global ncbi_go_file; global PathDir; global remove_download_files; remove_download_files = 'no'
    global system_list; system_list=[]; global system_codes; system_codes={}; global option_db; global option_group_db
    global file_location_defaults; global integrate_online_species; integrate_online_species = 'no'; global PathFile

    na = 'NA'; log = 'log'; no = 'no';
    run_mappfinder=no; mod=na; permutation=0; modifyDBs=no; overwrite_affycsv = na; overwrite_entrezgo = na
    incorporate_previous_associations=no; process_go=no; z_threshold=1.96; p_val_threshold=0.05; resources_to_analyze = ''
    change_threshold=2; max_member_count = 10000; output_dir = ''; input_dir = ''; denom_dir = ''; custom_sets_dir = ''
    ORA_algorithm = no; run_from_scratch = None; returnPathways = None
    
    option_group_db,option_db = importUserOptions('options')  ##Initially used to just get the info for species and array_type
    import copy
    importSpeciesInfo()
    if len(species_codes) == 0:
        integrate_online_species = 'yes'
        importSpeciesInfo() ### Occurs when an unknown error erases the species file

    file_location_defaults = importDefaultFileLocations()
    ncbi_go_file = file_location_defaults['EntrezGO'].Location()
    null,system_list,mod_list = importSystemInfo()
    system_code_db = {}
    for sysname in system_codes: ad = system_codes[sysname]; system_code_db[ad.SystemCode()]=sysname   
    species_code_db = {}
    for species in species_codes: sd = species_codes[species]; species_code_db[sd.SpeciesCode()]=species
    
    try: PathDir = file_location_defaults['PathDir'].Location()
    except Exception: PathDir = ''
    try: PathFile = file_location_defaults['PathFile'].Location()
    except Exception: PathFile = ''

    elite_db_versions=[]
    try: elite_db_versions = returnDirectoriesNoReplace('/Databases')
    except Exception:
        try: elite_db_versions=[]; 
        except Exception: null=[]### directory already exists      
    try: gene_database_dir = unique.getCurrentGeneDatabaseVersion()
    except Exception: gene_database_dir=''
    
    global gu
    ###Update this informatin in option_db which will be over-written after the user selects a species and array_type                
    current_species_names = getSpeciesList()  ### returns a list of the installed species names (full)
    species_list_augmented = copy.deepcopy(current_species_names)
    species_list_augmented2 = copy.deepcopy(species_list)
    
    option_db['species'].setAnalysisOptions(current_species_names)
    option_db['mod'].setAnalysisOptions(mod_list)
    resource_list = importResourceList()
    species_list.sort()
                        
    if len(elite_db_versions)>1:
        option_db['dbase_version'].setAnalysisOptions(elite_db_versions)
        option_db['dbase_version'].setDefaultOption(gene_database_dir)
    else:
        ### Otherwise, remove this option
        del option_db['dbase_version']; option_group_db['ORA'] = option_group_db['ORA'][1:]

    species_list_augmented.sort(); species_list_augmented2.sort()
    if len(species_list_augmented2)>0 and species_list_augmented2 != ['']:
        species_list_augmented.append('----'); species_list_augmented2.append('----'); species_list_augmented.append('New Species'); species_list_augmented.append('all-supported')
    else: species_list_augmented=[]; species_list_augmented2=[]
    if len(species_list_augmented)==0: species_list_augmented.append('New Species')
    species_list_augmented2.append('New Species')
    species_list_augmented.reverse()
    option_db['species_resources'].setAnalysisOptions(current_species_names) #species_list_augmented2[:-2]
    option_db['species_affy_update'].setAnalysisOptions(species_list_augmented)
    option_db['species_eg_update'].setAnalysisOptions(species_list_augmented)
    option_db['species_update'].setAnalysisOptions(species_list_augmented2)
        
    if run_parameter[0] == 'ORA': ### Occurs when selecting "Back" from Elite parameter window
        old_options = run_parameter[1]
        run_mappfinder = 'yes'; run_parameter = 'ORA'
        for option in old_options: ### Set options to user selected
            option_db[option].setDefaultOption(old_options[option])
        
    try:
        if run_parameter == 'skip' or run_parameter == 'intro' or run_parameter == 'ORA':
            if run_parameter!= 'ORA':
                root = Tk()
                root.title('GO-Elite: Main Dataset Parameters')
                gu = GUI(root,option_db,option_group_db['parameters'])
                run_from_scratch = gu.Results()['run_from_scratch']    

                if run_from_scratch == 'Analyze ID Lists': run_mappfinder = 'yes'; update_dbs = 'no'
                if run_from_scratch == 'Prune Existing Results': run_mappfinder = 'no'; update_dbs = 'no'
                if run_from_scratch == 'Update or Add Databases': run_mappfinder = 'no'; update_dbs= 'yes'
                if run_from_scratch == 'View Data on Pathways': run_mappfinder = 'NA'; update_dbs = 'NA'
                if update_dbs == 'no':
                    if len(species_codes)==0 or gene_database_dir == '' or len(elite_db_versions)==0: ### No species stored
                        print_out = "No species databases found. Select\ncontinue to proceed with species download."
                        IndicatorWindow(print_out,'Continue')
                        integrate_online_species = 'yes'
                        addOnlineSpeciesDatabases()
                        GO_Elite.importGOEliteParameters(['ORA',{}]); sys.exit()

            if run_mappfinder == 'yes':                        
                proceed = 'no'; update_dbs = 'no'
                while proceed == 'no':
                    root = Tk(); 
                    root.title('GO-Elite: Over-representation Analysis (ORA) Parameters')
                    gu = GUI(root,option_db,option_group_db['ORA'])
                    species = gu.Results()['species']
                    species_code = species_codes[species].SpeciesCode()
                    permutation = gu.Results()['permutation']
                    mod = gu.Results()['mod']
                    try: input_dir = gu.Results()['input_dir']
                    except KeyError: input_dir = ''
                    try: denom_dir = gu.Results()['denom_dir']
                    except KeyError: denom_dir = ''
                    try: custom_sets_dir = gu.Results()['custom_sets_dir']
                    except KeyError: custom_sets_dir = ''
                    try: ORA_algorithm = gu.Results()['ORA_algorithm']
                    except KeyError: ORA_algorithm = ''
                    for option in gu.Results(): ### Set options to user selected if error occurs
                        option_db[option].setDefaultOption(gu.Results()[option])
                    try:
                        geneGO_import_dir = '/Databases/'+species_code+'/gene-go'
                        gg = gene_associations.GrabFiles(); gg.setdirectory(geneGO_import_dir)
                        filedir,goname = gg.searchdirectory(mod)
                    except Exception: goname = ''
                    try: 
                        geneMAPP_import_dir = '/Databases/'+species_code+'/gene-mapp'
                        gm = gene_associations.GrabFiles(); gm.setdirectory(geneMAPP_import_dir)
                        filedir,mappname = gm.searchdirectory(mod)
                    except Exception: mappname = ''
                    if goname=='' and mappname=='':
                        print_out = "The primary system (aka MOD) %s is currently unavailable for that species." % mod
                        IndicatorWindow(print_out,'Continue')
                    elif len(input_dir)>0 and len(denom_dir)>0:
                        try:
                            null = int(permutation)
                            proceed = 'yes'
                        except Exception:
                            print_out = "Invalid numerical permutation entry. Try again."
                            IndicatorWindow(print_out,'Continue')
                        input_text_files = readDirText(input_dir)
                        input_denom_files = readDirText(denom_dir)
                        if len(input_text_files)>0 and len(input_denom_files)>0: proceed = 'yes'
                        else:
                            proceed = 'no'
                            if len(input_text_files)==0:
                                print_out = 'No files with the extension ".txt" found in the input directory.'
                                IndicatorWindow(print_out,'Continue')
                            else:
                                print_out = 'No files with the extension ".txt" found in the denominator directory.'
                                IndicatorWindow(print_out,'Continue')
                    else:
                        print_out = "Please designate an input and denominator file directory"
                        IndicatorWindow(print_out,'Continue')
        else:
            update_dbs = 'yes' ### Happens when a previous option was selected (or warning) that re-runs this interface to get a specific menu                
            modifyDBs = run_parameter

        if run_from_scratch == 'View Data on Pathways':
                root = Tk()
                root.title('GO-Elite: Visualize Data on WikiPathways')
                gu = GUI(root,'ViewWikiPathways',[])
        if update_dbs == 'yes':
            if run_parameter == 'skip' or run_parameter == 'intro':
                ###Update this informatin in option_db which will be over-written after the user selects a species and array_type
                root = Tk()
                root.title('GO-Elite: Update Options')
                #print option_list[i:i+1];kill
                gu = GUI(root,option_db,option_group_db['update'])
                modifyDBs = gu.Results()['modifyDBs1']
            if modifyDBs == 'Create/Modify Databases':
                root = Tk(); root.title('GO-Elite: Create/Modify/Update Databases')
                gu = GUI(root,option_db,option_group_db['update2'])
                modifyDBs = gu.Results()['modifyDBs2']

            if modifyDBs == 'Your Own Text Files':
                species = '---'
                while '--' in species:
                    root = Tk(); root.title('GO-Elite: Manually Add New Relationships')
                    gu = GUI(root,option_db,option_group_db['customOptions'])
                    species = gu.Results()['species_update'] 
                if species == 'New Species':
                    new_run_parameter = modifyDBs
                    modifyDBs = 'Add New Species Support'
                    update_options = ''
                else:
                    species_code = species_codes[species].SpeciesCode()
                    update_options = gu.Results()['update_options']
                    current_species_dirs = unique.read_directory('/Databases')
                    if species_code not in current_species_dirs:
                        print_out = "Support for this species not downloaded yet. Select\ncontinue to proceed with species download."
                        IndicatorWindow(print_out,'Continue')
                        integrate_online_species = 'yes'
                        addOnlineSpeciesDatabases()
                        GO_Elite.importGOEliteParameters(['modifyDBs2',{}]); sys.exit()
                    if update_options == 'Add New Gene System':
                        new_system_name= ''; new_system_code = ''
                        while new_system_name == '' or new_system_code =='':
                            root = Tk(); root.title('GO-Elite: Add New Gene System')
                            gu = GUI(root,option_db,option_group_db['addCustomSystem'])
                            new_system_name = gu.Results()['new_system_name']
                            new_system_code = gu.Results()['new_system_code']
                            new_mod = gu.Results()['new_mod']
                            try: gene_system_file = gu.Results()['gene_system_file']
                            except Exception: gene_system_file = ''
                            if new_system_code in system_code_db:
                                if new_system_name == system_code_db[new_system_code]:
                                    print_out = 'The system code and name entered already exist.\nSelect replace to change the "MOD" status.'
                                    fw = FeedbackWindow(print_out,'Replace',"Don't Replace")
                                    choice = fw.ButtonSelection()['button']
                                    if choice != 'Replace': GO_Elite.importGOEliteParameters(modifyDBs); sys.exit()  
                                else:
                                    print_out = 'The system code entered already exists.\nSelect replace to system information.'
                                    fw = FeedbackWindow(print_out,'Replace',"Don't Replace")
                                    choice = fw.ButtonSelection()['button']
                                    if choice != 'Replace': GO_Elite.importGOEliteParameters(modifyDBs); sys.exit()
                            elif new_system_name in system_codes:
                                print_out = 'The system name entered already exists.\nSelect replace to system information.'
                                fw = FeedbackWindow(print_out,'Replace',"Don't Replace")
                                choice = fw.ButtonSelection()['button']
                                if choice != 'Replace': GO_Elite.importGOEliteParameters(modifyDBs); sys.exit()
                            if new_mod == 'yes':
                                new_mod = 'MOD'
                                if len(gene_system_file) == 0:
                                    new_system_name = ''; new_system_code = ''
                            else:
                                new_mod = ''
                            if '-' in new_system_name:
                                print_out = 'The system name contains the invalid character "-".\nPlease change and try again.'
                                fw = FeedbackWindow(print_out,'Continue',"")
                                choice = fw.ButtonSelection()['button']
                                new_system_name=''; new_system_code=''
                            elif new_system_name == '' or new_system_code =='': 
                                print_out = "No system name or system code specified."
                                if new_mod == 'MOD':
                                    print_out = "An ID annotation file for this system must be included (required for MODs)"
                                IndicatorWindow(print_out,'Continue')
                            
                        ### Add new system info
                        ad = SystemData(new_system_code,new_system_name,new_mod)
                        system_codes[new_system_name] = ad
                        exportSystemInfo()
                        if len(gene_system_file)>0:
                            root = Tk(); analysis = 'updateAnnotationFiles'
                            overwrite_annotations = 'yes'
                            values = gene_system_file,new_system_name,overwrite_annotations,species_code
                            StatusWindow(root,analysis,values)
                        GO_Elite.importGOEliteParameters(modifyDBs); sys.exit()                    
                    elif update_options == 'Add New Relationship Table':
                        proceed = 'no'
                        system_list.sort(); system_list.reverse()
                        system_list += ['----','Gene Ontology']+ resource_list[1:] + ['New Resource'] ### Add all existing resources to this list
                        system_list.reverse()
                        option_db['mod_to_update'].setAnalysisOptions(mod_list)
                        option_db['system_to_update'].setAnalysisOptions(system_list)                        
                        while proceed == 'no':
                            root = Tk(); root.title('GO-Elite: Add New Relationship Table')
                            gu = GUI(root,option_db,option_group_db['customRelationships'])
                            mod_to_update = gu.Results()['mod_to_update']
                            system_to_update = gu.Results()['system_to_update']
                            new_resource_name = gu.Results()['new_resource_name']
                            new_resource_type = gu.Results()['new_resource_type']
                            try: update_custom_relationship_file = gu.Results()['update_custom_relationship_file']
                            except Exception: update_custom_relationship_file = ''
                            overwrite_relationships = 'yes'
                            relationship_to_update = mod_to_update+'-'+system_to_update
                            if 'Gene Ontology' in system_to_update:
                                relationship_to_update = mod_to_update+'-GeneOntology'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-go/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'Disease Ontology' in system_to_update:
                                relationship_to_update = mod_to_update+'-CTDOntology'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-go/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'Phenotype Ontology' in system_to_update:
                                relationship_to_update = mod_to_update+'-MPhenoOntology'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-go/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'GOSlim' in system_to_update:
                                relationship_to_update = mod_to_update+'-GOSlim'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-go/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'WikiPathways' in system_to_update:
                                relationship_to_update = mod_to_update+'-MAPP'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-mapp/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'miRNA Targets' in system_to_update:
                                relationship_to_update = mod_to_update+'-microRNATargets'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-mapp/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'BioMarkers' in system_to_update:
                                relationship_to_update = mod_to_update+'-BioMarkers'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-mapp/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'Domains' in system_to_update:
                                relationship_to_update = mod_to_update+'-Domains'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-mapp/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'PathwayCommons' in system_to_update:
                                relationship_to_update = mod_to_update+'-PathwayCommons'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-mapp/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'KEGG' in system_to_update:
                                relationship_to_update = mod_to_update+'-KEGG'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-mapp/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'Transcription Factor Targets' in system_to_update:
                                relationship_to_update = mod_to_update+'-TFTargets'
                                relationship_to_update_file ='Databases/'+species_code+'/gene-mapp/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            elif 'New Resource' in system_to_update:
                                relationship_to_update = mod_to_update+'-'+new_resource_name
                                if new_resource_type == 'Ontology':
                                    relationship_to_update_file ='Databases/'+species_code+'/gene-go/'+relationship_to_update+'.txt'
                                else:
                                    relationship_to_update_file ='Databases/'+species_code+'/gene-mapp/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            else:
                                relationship_to_update_file ='Databases/'+species_code+'/uid-gene/'+relationship_to_update+'.txt'
                                file_present = verifyFile(relationship_to_update_file)
                            if file_present == 'yes':
                                print_out = 'The relationship table already exists. Only\nuse this menu to add new tables and use the\n"Update Existing Relationship Table" to update\nexisting relationships.'
                                IndicatorWindow(print_out,'Continue'); GO_Elite.importGOEliteParameters(modifyDBs); sys.exit()         
                            elif (mod_to_update != system_to_update) and ('----' not in system_to_update):
                                root = Tk(); analysis = 'updateRelationshipFiles'
                                values = update_custom_relationship_file,relationship_to_update_file,overwrite_relationships,species_code
                                StatusWindow(root,analysis,values)
                                proceed = 'yes'
                        GO_Elite.importGOEliteParameters(modifyDBs); sys.exit()
                    elif update_options == 'Update Existing Relationship Table':
                        import_dir1 = '/Databases/'+species_code+'/uid-gene'
                        import_dir2 = '/Databases/'+species_code+'/gene-mapp'
                        import_dir3 = '/Databases/'+species_code+'/gene-go'
                        #import_dir4 = '/Databases/'+species_code+'/gene'
                        proceed = 'yes'; systems_not_present = []
                        try: uid_gene_list = read_directory(import_dir1); uid_gene_list.sort()
                        except Exception: proceed = 'no'; systems_not_present.append('uid-gene'); uid_gene_list = []
                        try: gene_mapp_list = read_directory(import_dir2); gene_mapp_list.sort()
                        except Exception: systems_not_present.append('gene-mapp'); gene_mapp_list=[]
                        try: gene_go_list = read_directory(import_dir3); gene_go_list.sort()
                        except Exception: systems_not_present.append('gene-go'); gene_go_list=[]
                        #try: gene_list = read_directory(import_dir4); gene_list.sort()
                        #except Exception: proceed = 'no'; systems_not_present.append('gene')
                        if len(gene_go_list) == 1 and 'Ensembl_version.txt' in gene_go_list: gene_go_list=[]
                        #if len(gene_list) == 0: proceed = 'no'; systems_not_present.append('gene')
                        elif len(uid_gene_list) == 0: proceed = 'no'; systems_not_present.append('uid-gene')
                        elif len(gene_mapp_list) == 0 and len(gene_go_list) == 0: proceed = 'no'; systems_not_present.append('gene-mapp'); systems_not_present.append('gene-go')
                        if (proceed == 'no') or ('gene-mapp' in systems_not_present and 'gene-go' in systems_not_present):
                            print_out = "Please note: The species directory does not appear\nto have a valid %s table(s). You may need to add before proceeding." % str(systems_not_present)[1:-1]
                            IndicatorWindow(print_out,'Continue')
                            #GO_Elite.importGOEliteParameters('skip'); sys.exit()
                        relationship_dir = gene_go_list+['----']+gene_mapp_list+['----']+uid_gene_list#+['----']+gene_list
                        relationship_files=[]
                        for i in relationship_dir:
                            if '-' in i: relationship_files.append(string.replace(i,'.txt',''))
                        option_db['relationship_to_update'].setAnalysisOptions(relationship_files)
                        root = Tk(); root.title('GO-Elite: Update Existing Relationship Table')
                        gu = GUI(root,option_db,option_group_db['customUpdate'])
                        relationship_to_update = gu.Results()['relationship_to_update']
                        update_relationship_file = gu.Results()['update_relationship_file']
                        overwrite_relationships = gu.Results()['overwrite_relationships']
                        
                        if relationship_to_update != '----':
                            file = relationship_to_update+'.txt'
                            if file in gene_mapp_list:
                                relationship_to_update = 'Databases/'+species_code+'/gene-mapp/'+relationship_to_update+'.txt'
                            elif file in gene_go_list:
                                relationship_to_update = 'Databases/'+species_code+'/gene-go/'+relationship_to_update+'.txt'
                            root = Tk(); analysis = 'updateRelationshipFiles'
                            values = update_relationship_file,relationship_to_update,overwrite_relationships,species_code
                            StatusWindow(root,analysis,values)
            elif modifyDBs == 'Manually Add New Relationships':
                print_out = "Please Install a Species Database or Add New Species Support"
                IndicatorWindow(print_out,'Continue')
                GO_Elite.importGOEliteParameters('skip'); sys.exit()  
            if modifyDBs == 'EntrezGene-GO Associations':
                if len(species_codes)>0:
                    species = '---'
                    while '--' in species:
                        root = Tk(); root.title('GO-Elite: EntrezGene-GO Update Options')
                        gu = GUI(root,option_db,option_group_db['updateEntrezGO'])
                        species = gu.Results()['species_eg_update']; compatible_mods = []
                    try: species_code = [species_codes[species].SpeciesCode()]
                    except KeyError: species_code = '--'
                    download_entrez_go = gu.Results()['download_entrez_go']
                    incorporate_previous_associations_EG = gu.Results()['incorporate_previous_associations_EG']
                    if incorporate_previous_associations_EG == 'update previous relationships':
                        rewrite_existing_EG = 'no'
                    else: rewrite_existing_EG = 'yes'
                    #overwrite_entrezgo = gu.Results()['overwrite_entrezgo']
                    overwrite_entrezgo = 'over-write previous databases' ### Don't include this option... save to NewDatabases
                    remove_download_files = gu.Results()['delete_entrezgo']
                    status = ''
                    
                    if species_code != '----' and species != 'New Species' and species != 'all-supported':
                        taxid = species_codes[species].TaxID()
                        #species = 'New Species' ### Must designate a taxid
                        compatible_mods = species_codes[species].Systems()
                    if species == 'New Species':
                        new_run_parameter = modifyDBs
                        modifyDBs = 'Add New Species Support'; species=''
                    if species == 'all-supported':
                        tax_ids=[]; species_code = []
                        for species in species_codes:
                            taxid = species_codes[species].TaxID(); tax_ids.append(taxid)
                            sp_code = species_codes[species].SpeciesCode(); species_code.append(sp_code)
                    elif len(species)>1: tax_ids = [species_codes[species].TaxID()]
                    if modifyDBs != 'Add New Species Support':
                        if download_entrez_go == 'yes':
                            root = Tk(); analysis = 'DownloadEntrezGO'
                            values = ncbi_go_file,'BuildDBs/Entrez/Gene2GO/','txt',tax_ids,species_code,overwrite_entrezgo, rewrite_existing_EG
                            StatusWindow(root,analysis,values)
                        root = Tk(); analysis = 'EntrezGOExport'
                        values = tax_ids,species_code,status,option_db,option_group_db['gene2go_not_found'],overwrite_entrezgo, rewrite_existing_EG
                        StatusWindow(root,analysis,values)
                else:
                    print_out = "Please Install a Species Database or or Add New Species Support"
                    IndicatorWindow(print_out,'Continue')
                    GO_Elite.importGOEliteParameters(modifyDBs); sys.exit()  
            if modifyDBs == 'Download Species Databases':
                integrate_online_species = 'yes'
                addOnlineSpeciesDatabases()
                GO_Elite.importGOEliteParameters('skip'); sys.exit()
            if modifyDBs == 'Ensembl Associations':
                print_out = "This menu is intended to download and integrate Ensembl\nrelationships from versions and species NOT supported\nby GO-Elite. Current GO-Elite databases already contains\nall relationships for supported species and versions."
                IndicatorWindow(print_out,'Continue')
                
                global process_Ens_go
                ### Download Ensembl Species Information from FTP Server
                root = Tk(); analysis = 'getCurrentEnsemblSpecies'
                values = None
                StatusWindow(root,analysis,values)
                option_db['species_ensembl_update'].setAnalysisOptions(ensembl_species)
                option_db['ensembl_version'].setAnalysisOptions(ensembl_versions)
                root = Tk(); root.title('GO-Elite: Extract Relationships from Online Ensembl Databases')
                gu = GUI(root,option_db,option_group_db['updateEnsemblSQL'])
                species = gu.Results()['species_ensembl_update']
                ensembl_version = gu.Results()['ensembl_version']
                try: species_code = species_codes[species].SpeciesCode()
                except KeyError: species_code = '--'
                incorporate_previous_ensembl_associations = gu.Results()['incorporate_previous_ensembl_associations']
                if incorporate_previous_ensembl_associations == 'update previous relationships': rewrite_existing = 'no'
                else: rewrite_existing = 'yes'
                #overwrite_ensembl = gu.Results()['overwrite_ensembl']
                overwrite_ensembl = 'over-write previous databases' ### Don't include this option... save to NewDatabases
                remove_download_files = gu.Results()['delete_ensembl']
                force = gu.Results()['download_latest_ensembl_files']
                ### Add species information to database
                if species_code == '--':
                    genus,species_code = string.split(species,' '); species_code = genus[0]+species_code[0]
                    species_taxid = '' ### Ideally should have but can assign when updating EntrezGene directly
                    compatible_mods = ['En']
                    sd = SpeciesData(species_code,species,compatible_mods,species_taxid)
                    species_codes[species] = sd
                    exportSpeciesInfo(species_codes)

                root = Tk(); analysis = 'getVersionedEnsExternalDB'
                values = species,ensembl_version
                StatusWindow(root,analysis,values)

                external_dbs, external_system, array_systems, external_ids = importExternalDBs(species)
                option_db['include_ens1'].setAnalysisOptions(external_dbs)
                option_db['include_ens2'].setAnalysisOptions(external_dbs)
                option_db['include_ens3'].setAnalysisOptions(external_dbs)
                option_db['include_ens4'].setAnalysisOptions(external_dbs)
                option_db['include_ens5'].setAnalysisOptions(external_dbs)
                root = Tk(); root.title('GO-Elite: Extract Relationships from Online Ensembl Databases')
                gu = GUI(root,option_db,option_group_db['updateEnsemblSQL2'])
                include_ens1 = gu.Results()['include_ens1']
                include_ens2 = gu.Results()['include_ens2']
                include_ens3 = gu.Results()['include_ens3']
                include_ens4 = gu.Results()['include_ens4']
                include_ens5 = gu.Results()['include_ens5']
                process_Ens_go = gu.Results()['process_Ens_go']
                externalDBName_list = [include_ens1,include_ens2,include_ens3,include_ens4,include_ens5]
                if species_code == 'Hs' or species_code == 'Mm' or species_code == 'Rn':
                    for id in externalDBName_list:
                        if id in array_db:
                            print_out = "Attention: Integrating array information for\nthis species will increase download time and will\nrequire > 2GB of hard-drive space. Please confirm\nbefore continuing."
                            IndicatorWindow(print_out,'Continue')
                
                root = Tk(); analysis = 'EnsemblSQLImport'
                values = species_code,species,child_dirs,externalDBName_list,overwrite_ensembl,rewrite_existing,force,ensembl_version,external_system
  
                for external_db in externalDBName_list:
                    if len(external_db)>1:
                        syscode = external_system[external_db]
                        if syscode not in system_code_db:
                            ad = SystemData(syscode,external_db,'')
                            system_codes[external_db] = ad
                            system_code_db[ad.SystemCode()]= external_db
                exportSystemInfo()
                StatusWindow(root,analysis,values)
                
            if 'WikiPathways' in modifyDBs or 'Affymetrix Annotation files' in modifyDBs:
                if 'Affymetrix Annotation files' in modifyDBs:
                    print_out = "This menu is intended to integrate the most recent\n"
                    print_out+= "downloaded Affymetrix CSV annotation file relationships.\n"
                    print_out+= "However, all supported Affymetrix species CSV annotation\n"
                    print_out+= "   files relationships are present in the downloadable databases   \n"
                    print_out+= "from GO-Elite (e.g., EnsMart56Plus - where Plus indicates\n"
                    print_out+= "the addition of Affymetrix CSV relationships in addition to\nthose from Ensembl)."
                    IndicatorWindow(print_out,'Continue')
                try: option_db['species_affy_update'].setAnalysisOptions(species_list_augmented)
                except IOError: null = []
                if 'WikiPathways' in modifyDBs: title = 'GO-Elite: WikiPathways Update'
                else: title = 'GO-Elite: Affymetrix Annotation Update'
                
                if modifyDBs == 'Affymetrix Annotation files': analysis_option = 'updateAffyCSV'
                else: analysis_option = 'WikiPathways'
                species = '---'
                while '--' in species:
                    root = Tk(); root.title(title)
                    gu = GUI(root,option_db,option_group_db[analysis_option])
                    species = gu.Results()['species_affy_update']
                if modifyDBs == 'Affymetrix Annotation files':
                    process_go = gu.Results()['process_go']
                    parse_genesets = 'no' #gu.Results()['parse_genesets']
                    integrate_affy_associations = 'yes'
                #else: process_go = 'no'; parse_genesets = 'yes'; integrate_affy_associations = 'no'
                #overwrite_affycsv = gu.Results()['overwrite_affycsv']
                overwrite_affycsv = 'over-write previous databases' ### Don't include this option... save to NewDatabases
                try: incorporate_previous_associations = gu.Results()['incorporate_previous_associations']
                except KeyError: incorporate_previous_associations = gu.Results()['incorporate_previous_associations_WP']
                try: get_mapped_file_version = gu.Results()['get_mapped_file_version']
                except KeyError: get_mapped_file_version = ''
                if incorporate_previous_associations == 'update previous relationships':
                    incorporate_previous_associations = 'yes'
                else: incorporate_previous_associations = 'no'
                if species == 'all-supported':
                    species=[]; species_code = []
                    for sp in species_codes:
                        species.append(sp)
                        sp_code = species_codes[sp].SpeciesCode(); species_code.append(sp_code)
                elif species == 'New Species':
                    new_run_parameter = modifyDBs
                    modifyDBs = 'Add New Species Support'
                elif species == '----':
                    print_out = 'Please select a valid species.'
                    IndicatorWindow(print_out,'Continue')
                    GO_Elite.importGOEliteParameters(modifyDBs); sys.exit()  
                else:
                    species_code = [species_codes[species].SpeciesCode()]
                    species = [species]
                
                annotation_files_present = 'no'
                if modifyDBs != 'Add New Species Support':
                    while annotation_files_present == 'no':
                        try: import_dir = '/BuildDBs/Affymetrix/'+species_code[0]
                        except Exception:
                            print_out = "Please Install a Species Database or Aor Add New Species Support"
                            IndicatorWindow(print_out,'Continue')
                            GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()  
                        try: dir_list = read_directory(import_dir)  #send a sub_directory to a function to identify all files in a directory
                        except Exception: ### Occurs if species dir not found
                            dir_list = []
                        for affy_data in dir_list:    #loop through each file in the directory to output results
                            affy_data_dir = import_dir[1:]+'/'+affy_data
                            if '.csv' in affy_data_dir: annotation_files_present = 'yes'
                        #meta_file_present = verifyFile('Databases/'+species[0]+'/'+'uid-gene/Ensembl_EntrezGene-meta.txt')
                        if modifyDBs == 'WikiPathways': annotation_files_present = 'yes' # and meta_file_present == 'yes'
                        elif annotation_files_present == 'no':
                            print_out = "No Affymetrix annotation files found in the directory:\n"+osfilepath(import_dir[1:])
                            print_out += "\n\nAtleast one CSV file is needed to parse Affymetrix gene relationships and subsequently"
                            print_out += "\nprocess WikiPathway relationships.\n\nPress Continue to select a folder containing CSV file(s) from your computer."
                            IndicatorWindow(print_out,'Continue'); assinged = 'no'
                            while assinged == 'no':
                                root = Tk()
                                root.title('GO-Elite: Select Affymetrix Annotation File(s)')
                                gu = GUI(root,option_db,option_group_db['InputAnnotationFolder'])
                                input_annotation_dir = gu.Results()['input_annotation_dir'] + '/'
                                dir_list = read_directory(input_annotation_dir)
                                for input_annotation_file in dir_list:
                                    input_annotation_file = input_annotation_dir+input_annotation_file
                                    input_annotation_lower = string.lower(input_annotation_file)
                                    if '.csv' in input_annotation_lower:
                                        assinged = 'yes'
                                        ###Thus the CSV was confirmed, so copy it over to BuildDB
                                        icf_list = string.split(input_annotation_file,'/'); csv_short = icf_list[-1]
                                        destination_parent = import_dir[1:]+'/'
                                        info_list = input_annotation_file,osfilepath(destination_parent+csv_short)
                                        StatusWindow(Tk(),'copy',info_list)
                    root = Tk(); analysis = 'AffyCSV'
                    values = system_codes,species_code,species,incorporate_previous_associations,process_go,parse_genesets,integrate_affy_associations,overwrite_affycsv,get_mapped_file_version
                    StatusWindow(root,analysis,values)
                    
            if modifyDBs == 'Add New Species Support':
                root = Tk(); root.title('GO-Elite: Add New Species Support')
                gu = GUI(root,option_db,option_group_db['NewSpecies'])
                species_code = gu.Results()['new_species_code']
                try: dbase_version = gu.Results()['dbase_version']; exportDBversion(db_version)
                except Exception: null = []
                new_species_name = gu.Results()['new_species_name']
                species_taxid = gu.Results()['species_taxid']
                compatible_mods =[]
                if species_code in species_code_db:
                    if new_species_name == species_code_db[species_code]:
                        print_out = 'The species code and name entered already exist.\nSelect replace to add/change the taxid.'
                        fw = FeedbackWindow(print_out,'Replace',"Don't Replace")
                        choice = fw.ButtonSelection()['button']
                        if choice != 'Replace': GO_Elite.importGOEliteParameters(new_run_parameter); sys.exit()  
                    else:
                        print_out = 'The species code entered already exists.\nSelect replace to change system information.'
                        fw = FeedbackWindow(print_out,'Replace',"Don't Replace")
                        choice = fw.ButtonSelection()['button']
                        if choice != 'Replace': GO_Elite.importGOEliteParameters(new_run_parameter); sys.exit()
                elif new_species_name in species_codes:
                    print_out = 'The species name entered already exists.\nSelect replace to change system information.'
                    fw = FeedbackWindow(print_out,'Replace',"Don't Replace")
                    choice = fw.ButtonSelection()['button']
                    if choice != 'Replace': GO_Elite.importGOEliteParameters(new_run_parameter); sys.exit()
                            
                sd = SpeciesData(species_code,new_species_name,compatible_mods,species_taxid)
                species_codes[new_species_name] = sd
                exportSpeciesInfo(species_codes)
                current_species_dirs = unique.read_directory('/Databases')
                if len(current_species_dirs) == 0:
                    ### Hence no official database has been added yet and adding a novel species
                    new_dir = 'Databases/EnsMart00/'+species_code
                    export.createExportFolder(new_dir) ### Create a species directory
                    export.createExportFolder(new_dir+'/gene') ### Add default directories
                    export.createExportFolder(new_dir+'/gene-go') ### Add default directories
                    export.createExportFolder(new_dir+'/uid-gene') ### Add default directories
                    export.createExportFolder(new_dir+'/gene-mapp') ### Add default directories
                    exportDBversion('EnsMart00')
                else:
                    ### An EnsMart directory thus exists
                    new_dir = 'Databases/'+species_code
                    try: export.createExportFolder(new_dir) ###Re-Create directory if deleted
                    except Exception: null=[]
                ### Indicate that support for this species is now added
                print_out = new_species_name+" succesfully added to database. Proceed\nwith addition of ID systems and relationships."
                ProcessCompletedWindow(print_out,'Process Completed')
                GO_Elite.importGOEliteParameters(new_run_parameter); sys.exit()
                
            if modifyDBs == 'Ontology structure':
                print_out = 'This menu is intended to download Ontology Structure annotations\nthat are NOT SUPPORTED from the "Database Species Download" menu.\nHowever, if you would like to get the very most recent files, proceed.'
                IndicatorWindow(print_out,'Continue')
                root = Tk(); root.title('GO-Elite: Update Options')
                gu = GUI(root,option_db,option_group_db['updateOBO'])
                OBO_url = gu.Results()['OBO_url']
                
                update_OBO = 'yes' #gu.Results()['update_OBO']

                root = Tk(); analysis = 'UpdateOBO'
                values = file_location_defaults,update_OBO,OBO_url
                StatusWindow(root,analysis,values)
                
            if modifyDBs == 'Additional Resources':
                ### Allow update of any supported resource (e.g., Ontologies, WikiPathways)
                current_species_dirs = unique.read_directory('/Databases')
                
                if len(current_species_dirs) == 0 or 'EnsMart' in current_species_dirs[0]: ### Either no species present or problem with the current config
                    print_out = "No species support currently found. Select\ncontinue to proceed with species download."
                    IndicatorWindow(print_out,'Continue')
                    integrate_online_species = 'yes'
                    addOnlineSpeciesDatabases()
                    GO_Elite.importGOEliteParameters(['modifyDBs2',{}]); sys.exit()
                    
                option_db['resource_to_update'].setAnalysisOptions(resource_list)
                root = Tk(); root.title('GO-Elite: Resource Update Options')
                gu = GUI(root,option_db,option_group_db['AdditionalResources'])
                species = gu.Results()['species_resources']
                additional_resources = gu.Results()['resource_to_update']
                species_code = species_codes[species].SpeciesCode()
                current_species_dirs = unique.read_directory('/Databases')
                root = Tk(); analysis = 'getAdditionalOnlineResources'
                values = species_code,additional_resources
                StatusWindow(root,analysis,values)
                
        if update_dbs == 'no':
            proceed = 'no'
            while proceed == 'no':
                try: all_species_codes = [species_code]
                except Exception:
                    all_species_codes=[]
                    for species in current_species_names:
                        all_species_codes.append(species_codes[species].SpeciesCode())
                for species_code in all_species_codes:
                    ### Augment the default resources to filter with others present
                    default_resources = option_db['resources_to_analyze'].AnalysisOptions()
                    #print option_db['species'].setAnalysisOptions(current_species_names)
                    import_dir1 = '/Databases/'+species_code+'/gene-mapp'
                    import_dir2 = '/Databases/'+species_code+'/gene-go'
                    try:
                        gene_mapp_list = read_directory(import_dir1)
                        gene_mapp_list.sort()
                        for file in gene_mapp_list:
                            resource = string.split(file,'-')[-1][:-4]
                            if resource != 'MAPP' and resource not in default_resources and '.txt' in file:
                                default_resources.append(resource)
                    except Exception: null=[]
                    try:
                        gene_go_list = read_directory(import_dir2)
                        gene_go_list.sort()
                        for file in gene_go_list:
                            resource = string.split(file,'-')[-1][:-4]
                            if resource != 'GeneOntology' and resource not in default_resources and 'version' not in resource and '.txt' in file:
                                default_resources.append(resource)
                    except Exception: null=[]
                
                option_db['resources_to_analyze'].setAnalysisOptions(default_resources)
                root = Tk(); root.title('GO-Elite: Redundancy Filtering Options')
                gu = GUI(root,option_db,option_group_db['elite'])
                z_threshold = gu.Results()['z_threshold']
                resources_to_analyze = gu.Results()['resources_to_analyze']
                p_val_threshold = gu.Results()['p_val_threshold']
                #include_headers_in_output = gu.Results()['include_headers_in_output']
                filter_method = gu.Results()['filter_method']
                try:
                    returnPathways = gu.Results()['returnPathways']
                    if returnPathways == 'no': returnPathways = None
                except Exception: returnPathways = None
                try: output_dir = gu.Results()['output_dir']
                except KeyError: output_dir = ''
                try: input_dir = gu.Results()['mappfinder_dir']
                except KeyError: null = ''
                if len(input_dir)>0 and len(output_dir)==0:
                    print_out = "Please select an output directory."
                    IndicatorWindow(print_out,'Continue')
                else:
                    try:
                        change_threshold = int(gu.Results()['change_threshold'])-1
                        null = float(z_threshold)
                        null = float(p_val_threshold)
                        proceed = 'yes'
                    except Exception:
                        print_out = "Invalid numerical entry. Try again."
                        IndicatorWindow(print_out,'Continue')
                    try:
                        max_member_count = int(gu.Results()['max_member_count'])
                        if max_member_count < 1: max_member_count = 10000
                    except Exception: null=[]
    except AttributeError,e:
        print 'Uknown error encountered... GO-Elite is exiting'
        print_out = e; "Unknown User Interface Error Encoutered"
        WarningWindow(print_out,print_out); root.destroy(); sys.exit()
    file_dirs = input_dir, denom_dir, output_dir, custom_sets_dir
    if ORA_algorithm == 'Fisher Exact Test':
        permutation = 'FisherExactTest'
    else:
        try: permutation = int(permutation)
        except Exception: null=[]
    return species, run_mappfinder, mod, permutation, filter_method, z_threshold, p_val_threshold, change_threshold, resources_to_analyze, max_member_count, returnPathways, file_dirs

def verifyFile(filename):
    fn=filepath(filename); file_found = 'yes'
    try:
        for line in open(fn,'rU').xreadlines():break
    except Exception: file_found = 'no'
    return file_found

def importResourceList():
    filename = 'Config/resource_list.txt'
    fn=filepath(filename); resource_list=[]
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        resource = data
        resource_list.append(resource)
    return resource_list

def updateOBOfiles(file_location_defaults,update_OBO,OBO_url,root):
    run_parameter = "Create/Modify Databases"
    if update_OBO == 'yes':
        import OBO_import
        c = OBO_import.GrabFiles()
        c.setdirectory('/OBO'); file_dirs = c.searchdirectory('.ontology')+c.searchdirectory('.obo')
        if len(OBO_url)>0: obo = OBO_url
        else: ### If not present, get the gene-ontology default OBO file
            obo = file_location_defaults['OBO'].Location()
        fln,status = update.download(obo,'OBO/','')
        run_parameter='Create/Modify Databases'
        if 'Internet' not in status:
            OBO_import.moveOntologyToArchiveDir()

            print_out = 'Finished downloading the latest Ontology OBO files.'
            print print_out

            try: system_codes,source_types,mod_types = GO_Elite.getSourceData()
            except Exception: null=[]
            
            if root !='' and root !=None:
                InfoWindow(print_out,'Update Complete!')
                continue_to_next_win = Button(text = 'Continue', command = root.destroy)
                continue_to_next_win.pack(side = 'right', padx = 10, pady = 10); root.mainloop()            
                GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
            else: null=[]
        else:
            if root !='' and root !=None: WarningWindow(status,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
            else: print status
    else:
        print_out = 'Download Aborted.'
        if root !='' and root !=None: WarningWindow(print_out,print_out); root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
        else: print print_out

def getOnlineEliteConfig(file_location_defaults,root):
    base_url = file_location_defaults['url'].Location()
    #fln,status = update.download(base_url+'Databases/','Databases/','')
    fln1,status1 = update.download(base_url+'Config/species_all.txt','Config/','')
    try:
        if 'Internet' not in status1:
            fln2,status2 = update.download(base_url+'Config/source_data.txt','Config/','')
            fln3,status3 = update.download(base_url+'Config/versions.txt','Config/','')
            #fln4,status4 = update.download(base_url+'Config/resource_list.txt','Config/','')
            print 'Finished downloading the latest configuration files.'; root.destroy()
        else:
            run_parameter = 'skip'
            try: WarningWindow(status3,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
            except Exception: print status3; root.destroy(); sys.exit()
    except Exception: null=[]

def buildInferrenceTables(species_code):
    try: gene_associations.swapAndExportSystems(species_code,'Ensembl','EntrezGene') ### Allows for analysis of Ensembl IDs with EntrezGene based GO annotations (which can vary from Ensembl)
    except Exception: null=[] ### Occurs if EntrezGene not supported
    try: gene_associations.augmentEnsemblGO(species_code)
    except Exception: null=[] ### Occurs if EntrezGene not supported

    ### Build out these symbol association files
    try: gene_associations.importGeneData(species_code,('export','Ensembl'))
    except Exception: null=[] ### Occurs if EntrezGene not supported
    try: gene_associations.importGeneData(species_code,('export','EntrezGene'))
    except Exception: null=[] ### Occurs if EntrezGene not supported
    
def getAdditionalOnlineResources(species_code,additional_resources,root):
    if additional_resources == 'All Resources':
        additional_resources = importResourceList()
    else: additional_resources = [additional_resources]
    try:
        print 'Adding supplemental GeneSet and Ontology Collections'
        import GeneSetDownloader; force = 'yes'
        GeneSetDownloader.buildAccessoryPathwayDatabases([species_code],additional_resources,force)
        print_out = 'Finished incorporating additional resources.'
    except Exception:
        print_out = 'Download error encountered for additional ontologies and gene-sets...\nplease try again later.'
    InfoWindow(print_out,'Continue')
    try: root.destroy()
    except Exception: null=[]
    GO_Elite.importGOEliteParameters('skip'); sys.exit()
        
def getOnlineEliteDatabase(file_location_defaults,db_version,new_species_codes,download_obo,additional_resources,root):
    base_url = file_location_defaults['url'].Location()
    if download_obo == 'yes':
        fln,status = update.download(base_url+'Databases/'+db_version+'/OBO.zip','','')
        
    for species_code in new_species_codes:
        #print [base_url+'Databases/'+db_version+'/'+species_code+'.zip']
        fln,status = update.download(base_url+'Databases/'+db_version+'/'+species_code+'.zip','Databases/','')
        buildInferrenceTables(species_code)
        ### Attempt to download additional Ontologies and GeneSets
        if additional_resources != '---':
            if additional_resources == 'All Resources': additionalResources = importResourceList()
            else: additionalResources = [additional_resources]
            try:
                print 'Adding supplemental GeneSet and Ontology Collections'
                import GeneSetDownloader; force = 'yes'
                GeneSetDownloader.buildAccessoryPathwayDatabases([species_code],additionalResources,force)
            except Exception: print 'Download error encountered for additional ontologies and gene-sets...\nplease try again later.'
        try: buildNestedOntologies(species_code)
        except Exception: None
    if 'Internet' not in status:                
        if len(new_species_codes)>0:
            print 'Finished downloading',db_version,'species database files.'
            print_out = "New species data succesfully added to database."
        else:
            print 'Finished downloading',db_version,'OBO database files.'
            print_out = "Ontology structure data succesfully added to database."     
        InfoWindow(print_out,'Continue')            
        try: root.destroy()
        except Exception: null=[]
    else:
        if root !='' and root !=None: WarningWindow(status,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('skip'); sys.exit()
        else: print status; root.destroy(); sys.exit()

def buildNestedOntologies(species_code):
    current_species_dirs = unique.returnDirectories('/Databases')
    if species_code in current_species_dirs:
        try: export.deleteFolder(filepath('Databases/'+species_code+'/nested')) ### Delete the existing nested folder (will force rebuilding it)
        except Exception: null=[]
        ### Creates a nested GO (and stores objects in memory, but not needed
        system_codes,source_types,mod_types = GO_Elite.getSourceData()
        #print species_code,'Building Nested for mod types:',mod_types
        avaialble_ontologies = OBO_import.findAvailableOntologies(species_code,mod_types)
        for ontology_type in avaialble_ontologies:
            full_path_db,path_id_to_goid,null = OBO_import.buildNestedOntologyAssociations(species_code,mod_types,ontology_type)
                            
def exportEntrezGO(tax_ids,species_codees,status,option_db,option_list,overwrite_entrezgo,rewrite_existing_EG,root):
    if 'over-write previous' in overwrite_entrezgo: overwrite_entrezgo = 'over-write previous' ### The config name is longer
    if 'Internet' not in status:
        import BuildAffymetrixAssociations
        index = 0
        print "Begining to parse NCBI Gene Ontology annotations..."
        for species_code in species_codees:
            tax_id = tax_ids[index]
            if len(tax_id)>0: print "Looking for %s EntrezGene to Gene Ontology associations" % species_code
            try: run_status = BuildAffymetrixAssociations.parseGene2GO(tax_id,species_code,overwrite_entrezgo,rewrite_existing_EG)
            except Exception: run_status = 'no'
            index+=1
        if run_status == 'run':
            print_out = 'Finished building EntrezGene-GeneOntology associations files.'
            if remove_download_files == 'yes':
                status = export.deleteFolder('BuildDBs/Entrez/Gene2GO'); print status
            InfoWindow(print_out,'Update Complete!')
            continue_to_next_win = Button(text = 'Continue', command = root.destroy)
            continue_to_next_win.pack(side = 'right', padx = 10, pady = 10); root.mainloop()
            GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
        else:
            print_out = 'Gene2GO file not found. Select download to obtain database prior to extraction'
            WarningWindow(print_out,'File Not Found!')            
            print '\nThe file "gene2go.txt" was not found.\n';
            root.destroy() ###Have to destroy status interface
            root = Tk(); root.title('GO-Elite: Update Options')
            gu = GUI(root,option_db,option_list)
            update_go_entrez = gu.Results()['update_go_entrez']
            if update_go_entrez == 'yes':
                root = Tk(); analysis = 'DownloadEntrezGO'
                values = ncbi_go_file,'BuildDBs/Entrez/Gene2GO/','txt',tax_ids,species_codees,overwrite_entrezgo,rewrite_existing_EG
                StatusWindow(root,'DownloadEntrezGO',values)
                #fln,status = update.download(ncbi_go_file,'BuildDBs/Entrez/Gene2GO/','txt')
                #exportEntrezGO(tax_id,species_code,status,option_db,option_list,root)
            else:
                #print_out = 'Gene2GO file not found. Select download to obtain database prior to extraction'
                #WarningWindow(print_out,'File Not Found!');  root.destroy();
                GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
    else: WarningWindow(status,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()

class importEnsemblSQL:
    def __init__(self,species,species_full,child_dirs,externalDBName_list,overwrite_previous,rewrite_existing,force,ensembl_version,external_system,root):
        self._parent = root; import EnsemblSQL
        if process_Ens_go == 'yes': externalDBName_list.append('GO')
        configType = 'Basic'; iteration=0; proceed = 'yes'

        ### Export Ensembl version information
        try: current_dirversion = unique.getCurrentGeneDatabaseVersion()
        except Exception: current_dirversion = ''
        dirversion = string.replace(ensembl_version,'release-','EnsMart')
        if dirversion not in current_dirversion:
            ### For example: if EnsMart56Plus is current and EnsMart56 is the selected version, keep EnsMart56Plus as current
            exportDBversion(dirversion)

        ### Instead of getting the current version info, get the specific version
        if ensembl_version != 'current':
            try:
                child_dirs, ensembl_species, ensembl_versions = EnsemblSQL.getCurrentEnsemblSpecies(ensembl_version)
            except Exception,e:
                print_out = "A internet connection could not be established.\nPlease fix the problem before proceeding."
                WarningWindow(e,'Warning!'); root.destroy(); GO_Elite.importGOEliteParameters(run_parameter); sys.exit()
            if species_full not in ensembl_species:
                print_out = 'The selected species is unavailable for\nthe selected version of Ensembl.'
                WarningWindow(print_out,'Species Unavailable!')
                root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
        if proceed == 'yes':
            ensembl_sql_dir,ensembl_sql_description_dir = child_dirs[species_full]
            for externalDBName in externalDBName_list:
                if externalDBName != ' ':
                    if force == 'yes' and iteration == 1: force = 'no'
                    import EnsemblSQL; reload(EnsemblSQL)
                    
                    if force == 'yes':
                        output_dir = 'BuildDBs/EnsemblSQL/'+species+'/'
                        try:
                            if force == 'yes': ### Delete any existing data in the destination directory that can muck up tables from a new Ensembl build
                                export.deleteFolder(output_dir)
                        except Exception: null=[]
                        
                    if externalDBName in array_db:
                        analysisType = 'FuncGen'
                        try: EnsemblSQL.buildGOEliteDBs(species,ensembl_sql_dir,ensembl_sql_description_dir,externalDBName,configType,analysisType,overwrite_previous,rewrite_existing,external_system,force); iteration+=1              
                        except Exception, e:
                            print traceback.format_exc()
                            print_out = 'This version of Ensembl appears to have critical incompatibilites with GO-Elite.\nDownload the latest version of GO-Elite or contact the development team'
                            WarningWindow(print_out,'Species Unavailable!')
                            root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit() 
                    else:
                        analysisType = 'GeneAndExternal'
                        try: EnsemblSQL.buildGOEliteDBs(species,ensembl_sql_dir,ensembl_sql_description_dir,externalDBName,configType,analysisType,overwrite_previous,rewrite_existing,external_system,force); iteration+=1
                        except Exception, e:
                            print_out = 'This version of Ensembl appears to have critical incompatibilites with GO-Elite.\nDownload the latest version of GO-Elite or contact the development team'
                            WarningWindow(print_out,'Species Unavailable!')
                            root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit() 
            if remove_download_files == 'yes': export.deleteFolder('BuildDBs/EnsemblSQL/'+species)
            print_out = 'Finished installing the selected Ensembl databases.'
            InfoWindow(print_out,'Update Complete!')
        continue_to_next_win = Button(text = 'Continue', command = root.destroy)
        continue_to_next_win.pack(side = 'right', padx = 10, pady = 10);
        quit_button = Button(root,text='Quit', command=self.quit)
        quit_button.pack(side = 'right', padx = 10, pady = 10)
        root.mainloop()            
        GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
        
    def quit(self):
        #print "quit starts"
        #print "cleaning up things..."
        self._parent.quit()
        self._parent.destroy()
        sys.exit()

def importExternalDBs(species_full):
    filename = 'Config/EnsExternalDBs.txt'
    fn=filepath(filename); x = 0; external_dbs=[]; external_system={}; all_databases={}; external_ids={}
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        if x==0: x=1
        else:
            id, database, species_specific, exclude, system_code = string.split(data,'\t')
            external_ids[database] = int(id)
            if database != 'GO':
                all_databases[database]=system_code
                if (species_full == species_specific) or len(species_specific)<2:
                    if len(exclude)<2:
                        external_system[database] = system_code

    filename = 'Config/external_db.txt'; external_system2={}
    fn=filepath(filename)
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        try:
            t = string.split(data,'\t'); id = int(t[0]); database = t[1]
            external_ids[database] = id
            if database in external_system:
                external_system2[database] = external_system[database]
            elif database not in all_databases: ### Add it if it's new
                try:
                    try: system = database[:3]
                    except Exception: system = database[:2]
                    external_system2[database] = system
                except Exception: null=[]
        except Exception: null=[] ### Occurs when a bad end of line is present

    filename = 'Config/array.txt'
    global array_db; array_db={}
    fn=filepath(filename)
    for line in open(fn,'rU').readlines():             
        data = cleanUpLine(line)
        t = string.split(data,'\t')
        try:
            array = t[1]; vendor = t[3]
            database = vendor+'_'+array; array_db[database]=[]
            if database in external_system:
                external_system2[database] = external_system[database]
            if database in all_databases:
                external_system2[database] = all_databases[database]
            elif database not in all_databases: ### Add it if it's new
                try:
                    if vendor == 'AFFY': system = 'X'
                    if vendor == 'ILLUMINA': system = 'Il'
                    if vendor == 'CODELINK': system = 'Co'
                    if vendor == 'AGILENT': system = 'Ag'
                    else: system = 'Ma'  ###Miscelaneous Array type
                    external_system2[database] = system
                except Exception: null=[]
        except Exception: null=[]
    external_system = external_system2
    #try: del external_system['GO']
    #except Exception: null=[]
    for database in external_system: external_dbs.append(database)
    external_dbs.append(' '); external_dbs = unique.unique(external_dbs); external_dbs.sort()
    return external_dbs,external_system,array_db,external_ids

def downloadFiles(file,dir,file_type,tax_id,species_code,overwrite_entrezgo,rewrite_existing_EG,root):
    fln,status = update.download(file,dir,file_type)
    if 'Internet' in status:
        WarningWindow(status,'Error Encountered!'); root.destroy(); GO_Elite.importGOEliteParameters('Create/Modify Databases'); sys.exit()
    else:
        exportEntrezGO(tax_id,species_code,status,{},{},overwrite_entrezgo,rewrite_existing_EG,root)

def testPNGView():
        tl = Toplevel()
        sf = PmwFreeze.ScrolledFrame(tl,
                labelpos = 'n', label_text = '',
                usehullsize = 1, hull_width = 800, hull_height = 550)
        sf.pack(padx = 0, pady = 0, fill = 'both', expand = 1)
        frame = sf.interior()

        filename = "/Users/nsalomonis/Desktop/code/AltAnalyze/datasets/3'Array/Kristina-Athro/DataPlots/WP1403-GE.Homozygot_vs_Control_AMPK signaling.png"
        tl.title(filename)
        img = ImageTk.PhotoImage(file=filename)
        can = Canvas(frame)
        can.pack(fill=BOTH, padx = 0, pady = 0)
        w = img.width()
        h = height=img.height()
        
        can.config(width=w, height=h)        
        can.create_image(2, 2, image=img, anchor=NW)
        tl.mainloop()
        
if __name__ == '__main__':
    user_variables=[]
    run_parameter = 'intro'
    vals = getUserParameters(run_parameter)
    print vals
