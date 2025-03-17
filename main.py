import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import os
import sys
from datetime import datetime
from Bio import Entrez
import numpy as np

class LeishDomainsApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("LeishDomains")
        
        # Global variables
        self.theInputFiles = []
        self.choice = None
        self.choice2 = None
        
        # Initialize UI
        self.setup_ui()
        
    def setup_ui(self):
        """Setup the main UI components"""
        # Create main menu
        self.create_menu()
        
        # Create main frames
        self.create_main_frames()
        
        # Create buttons and widgets
        self.create_widgets()
        
    def create_menu(self):
        """Create the main menu bar"""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        
        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="New Session", command=self.new_session)
        file_menu.add_command(label="Load Session", command=self.load_session)
        file_menu.add_command(label="Save Session", command=self.save_session)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        # Tools menu
        tools_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Tools", menu=tools_menu)
        tools_menu.add_command(label="Control Panel", command=self.show_control_panel)
        
        # Help menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="About", command=self.show_about)
        
    def create_main_frames(self):
        """Create the main application frames"""
        # Main content frame
        self.main_frame = ttk.Frame(self.root, padding="10")
        self.main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Text area frame
        self.text_frame = ttk.Frame(self.main_frame)
        self.text_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Buttons frame
        self.buttons_frame = ttk.Frame(self.main_frame)
        self.buttons_frame.grid(row=1, column=0, sticky=(tk.W, tk.E))
        
    def create_widgets(self):
        """Create and setup all widgets"""
        # Text area
        self.text_area = tk.Text(self.text_frame, width=80, height=20)
        self.text_area.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Scrollbar for text area
        scrollbar = ttk.Scrollbar(self.text_frame, orient=tk.VERTICAL, command=self.text_area.yview)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.text_area['yscrollcommand'] = scrollbar.set
        
        # Buttons
        ttk.Button(self.buttons_frame, text="New Session", command=self.new_session).grid(row=0, column=0, padx=5)
        ttk.Button(self.buttons_frame, text="Load Session", command=self.load_session).grid(row=0, column=1, padx=5)
        ttk.Button(self.buttons_frame, text="Save Session", command=self.save_session).grid(row=0, column=2, padx=5)
        
    def new_session(self):
        """Create a new session"""
        folder = filedialog.askdirectory(initialdir='.')
        if not folder:
            return
            
        print(f"FOLDER: {folder}")
        os.chdir(folder)
        fnames = glob("*.txt") + glob("*.xls")
        os.chdir(os.getcwd())
        
        for f in fnames:
            namebase = os.path.basename(f)
            namedir = folder.split("/")[-1]
            namewdir = namedir + "\\" + namebase
            self.theInputFiles.append(namewdir)
            
        messagebox.showinfo("New Session", f"Folder imported:\n{namewdir}")
        print(f"\nNew Session... {folder}\ncomplete!")
        
    def load_session(self):
        """Load an existing session"""
        print("\n##########    LOADING SESSION      ############")
        
        filename_input = filedialog.askopenfilename(
            title="Open File Session",
            filetypes=[("File TXT", ".txt"), ("All files", ".*")]
        )
        
        if not filename_input:
            messagebox.showinfo("Load session", "No Session Loaded!")
            return
            
        name = os.path.basename(filename_input)
        trig = name.split("_")
        
        if trig == "session" or "SESSION":
            with open(filename_input) as aFile:
                for line in aFile.readlines():
                    line = line.strip()
                    self.theInputFiles.append(line)
                    
            messagebox.showinfo("Load Session", f"Complete:\n{name}")
            print(f"\nSession Loading... {name}\nComplete!")
        else:
            messagebox.showinfo("Load Session", "No Session loaded...")
            
    def save_session(self):
        """Save the current session"""
        # Implementation from original save_session function
        pass
        
    def show_control_panel(self):
        """Show the control panel window"""
        # Implementation from original show_control_panel function
        pass
        
    def show_about(self):
        """Show the about dialog"""
        about_text = """LeishDomains
Version 1.55
A tool for analyzing protein domains in Leishmania species."""
        
        messagebox.showinfo("About", about_text)
        
    def run(self):
        """Start the application"""
        self.root.mainloop()

def main():
    app = LeishDomainsApp()
    app.run()

if __name__ == "__main__":
    main() 