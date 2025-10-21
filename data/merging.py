"""Data merging and joining operations."""

import os
import random
from collections import defaultdict
from typing import List, Dict, Tuple, Optional
from pathlib import Path


class DataMerger:
    """Handles merging and joining of data files."""
    
    def __init__(self):
        self.logger = None  # Will be set by caller if needed
    
    def _log(self, message: str):
        """Log message if logger is available."""
        if self.logger:
            self.logger.info(message)
        else:
            print(message)
    
    def _get_output_filename(self, file1: str, file2: str, operation: str, 
                            count1: int, count2: int) -> str:
        """Generate output filename based on input files and operation."""
        try:
            if "/" in file2:
                namef2 = file2.split("/")[-1]
            elif "\\" in file2:
                namef2 = file2.split("\\")[-1]
            else:
                namef2 = file2
            namef2 = namef2.replace('.txt', '')
        except:
            namef2 = file2.replace('.txt', '')
        
        namef1 = file2.replace('.txt', '')
        return f"{namef1}_{count1}_{namef2}_{operation}.txt"
    
    def merge_by_accession(self, file1: str, file2: str, output_path: Optional[str] = None) -> str:
        """Merge two files by accession number (first column)."""
        ac1_dict = defaultdict(list)
        ac2_dict = defaultdict(list)
        
        # Read file1
        with open(file1, 'r', encoding='utf-8') as f:
            header1 = f.readline().strip()
            for line in f:
                line = line.strip()
                if line:
                    items = line.split('\t')
                    if items:
                        ac1 = items[0]
                        ac1_dict[ac1].append(line)
        
        # Read file2
        with open(file2, 'r', encoding='utf-8') as f:
            header2 = f.readline().strip()
            for line in f:
                line = line.strip()
                if line:
                    items = line.split('\t')
                    if items:
                        ac2 = items[0]
                        ac2_dict[ac2].append(line)
        
        # Generate output filename
        if not output_path:
            output_path = self._get_output_filename(
                file1, file2, "merged", len(ac1_dict), len(ac2_dict)
            )
        
        # Write merged data
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(f"{header1}\t{header2}\n")
            
            count = 0
            pasted = []
            
            # Merge matching accessions
            for ac, lines1 in ac1_dict.items():
                if ac in ac2_dict:
                    for line1 in lines1:
                        for line2 in ac2_dict[ac]:
                            new_line = f"{line1}\t{line2}"
                            if new_line not in pasted:
                                pasted.append(new_line)
                                f.write(f"{new_line}\n")
                                count += 1
                else:
                    # Only in file1
                    for line1 in lines1:
                        if line1 not in pasted:
                            pasted.append(line1)
                            f.write(f"{line1}\n")
                            count += 1
            
            # Only in file2
            for ac, lines2 in ac2_dict.items():
                if ac not in ac1_dict:
                    for line2 in lines2:
                        if line2 not in pasted:
                            pasted.append(line2)
                            f.write(f"{line2}\n")
                            count += 1
        
        self._log(f"Merged {count} lines from {len(ac1_dict)} and {len(ac2_dict)} accessions")
        return output_path
    
    def match_by_accession(self, file1: str, file2: str, output_path: Optional[str] = None) -> str:
        """Match files by accession number, keeping only matching entries."""
        ac1_dict = defaultdict(list)
        ac2_dict = defaultdict(list)
        
        # Read file1
        with open(file1, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line:
                    items = line.split('\t')
                    if items:
                        ac1 = items[0]
                        ac1_dict[ac1].append(line)
        
        # Read file2
        with open(file2, 'r', encoding='utf-8') as f:
            header2 = f.readline().strip()
            for line in f:
                line = line.strip()
                if line:
                    items = line.split('\t')
                    if items:
                        ac2 = items[0]
                        ac2_dict[ac2].append(line)
        
        # Generate output filename
        if not output_path:
            output_path = self._get_output_filename(
                file1, file2, "match", len(ac1_dict), len(ac2_dict)
            )
        
        # Write matching data
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(f"{header2}\n")
            
            count = 0
            pasted = []
            
            for ac, lines1 in ac1_dict.items():
                if ac in ac2_dict:
                    for line1 in lines1:
                        if line1 not in pasted:
                            pasted.append(line1)
                            f.write(f"{line1}\n")
                            count += 1
        
        self._log(f"Matched {count} lines from {len(ac1_dict)} accessions")
        return output_path
    
    def join_lines(self, file1: str, file2: str, output_path: Optional[str] = None) -> str:
        """Join lines from two files, removing duplicates."""
        if not output_path:
            output_path = self._get_output_filename(file1, file2, "join", 0, 0)
        
        pasted = []
        count = 0
        
        with open(output_path, 'w', encoding='utf-8') as out:
            # Read file1
            with open(file1, 'r', encoding='utf-8') as f:
                header = f.readline().strip()
                out.write(f"{header}\n")
                
                for line in f:
                    line = line.strip()
                    if line and line not in pasted:
                        pasted.append(line)
                        out.write(f"{line}\n")
                        count += 1
            
            # Read file2
            with open(file2, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line and line not in pasted:
                        pasted.append(line)
                        out.write(f"{line}\n")
                        count += 1
        
        self._log(f"Joined {count} unique lines")
        return output_path
