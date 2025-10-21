import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import os
import sys
from datetime import datetime
from Bio import Entrez
import numpy as np
from glob import glob

from core.session import SessionManager
from core.logger import get_logger

class LeishDomainsApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("LeishDomains")
        
        # Global variables
        self.theInputFiles = []
        self.choice = None
        self.choice2 = None
        self.session = SessionManager()
        self.logger = get_logger()
        
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
            
        self.logger.info("Selected folder: %s", folder)
        os.chdir(folder)
        fnames = glob("*.txt") + glob("*.xls")
        os.chdir(os.getcwd())
        
        for f in fnames:
            namebase = os.path.basename(f)
            namedir = folder.split("/")[-1]
            namewdir = namedir + "\\" + namebase
            self.theInputFiles.append(namewdir)
            
        last = namewdir if fnames else folder
        messagebox.showinfo("New Session", f"Folder imported:\n{last}")
        self.logger.info("New session initialized with %d files", len(self.theInputFiles))
        
    def load_session(self):
        """Load an existing session"""
        self.logger.info("Loading session file")
        
        filename_input = filedialog.askopenfilename(
            title="Open File Session",
            filetypes=[("File TXT", ".txt"), ("All files", ".*")]
        )
        
        if not filename_input:
            messagebox.showinfo("Load session", "No Session Loaded!")
            return
            
        name = os.path.basename(filename_input)
        base_upper = name.upper()
        try:
            loaded = self.session.load(filename_input)
        except Exception as exc:
            self.logger.exception("Failed to load session: %s", exc)
            messagebox.showerror("Load Session", f"Failed to load: {name}")
            return

        if loaded:
            self.theInputFiles.extend(loaded)
            messagebox.showinfo("Load Session", f"Loaded {len(loaded)} entries from:\n{name}")
            self.logger.info("Session loaded: %s with %d entries", name, len(loaded))
        else:
            messagebox.showinfo("Load Session", "Session file was empty.")
            
    def save_session(self):
        """Save the current session"""
        if not self.theInputFiles:
            messagebox.showinfo("Save Session", "Nothing to save yet.")
            return

        target = filedialog.asksaveasfilename(
            title="Save Session As",
            defaultextension=".txt",
            filetypes=[("Text file", ".txt"), ("All files", ".*")]
        )
        if not target:
            return

        try:
            # If user picked a path, persist using our manager but preserve exact location
            directory = os.path.dirname(target)
            # Save to temp path in the chosen directory to get a standard name
            saved_path = self.session.save(self.theInputFiles, directory=directory)
            # If the chosen filename differs from default, move/overwrite
            if os.path.abspath(saved_path) != os.path.abspath(target):
                with open(saved_path, "r", encoding="utf-8") as src, open(target, "w", encoding="utf-8") as dst:
                    dst.write(src.read())
                try:
                    os.remove(saved_path)
                except Exception:
                    pass
            messagebox.showinfo("Save Session", f"Session saved:\n{target}")
            self.logger.info("Session saved: %s (%d entries)", target, len(self.theInputFiles))
        except Exception as exc:
            self.logger.exception("Failed to save session: %s", exc)
            messagebox.showerror("Save Session", "Failed to save session.")
        
    def show_control_panel(self):
        """Show the control panel window"""
        panel = tk.Toplevel(self.root)
        panel.title("Control Panel")
        ttk.Label(panel, text="Coming soon: analysis and tools").grid(row=0, column=0, padx=10, pady=10)
        
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