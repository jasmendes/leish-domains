"""Data filtering and search operations."""

import os
import random
from typing import List, Optional


class DataFilter:
    """Handles data filtering and search operations."""
    
    def __init__(self):
        self.logger = None  # Will be set by caller if needed
    
    def _log(self, message: str):
        """Log message if logger is available."""
        if self.logger:
            self.logger.info(message)
        else:
            print(message)
    
    def filter_by_uniprot_ids(self, uniprot_ids: List[str], filename: str, 
                             output_path: Optional[str] = None) -> str:
        """Filter file by list of UniProt IDs."""
        if not uniprot_ids:
            raise ValueError("No UniProt IDs provided")
        
        if not output_path:
            name = random.choice(uniprot_ids)
            output_path = f"{filename.replace('.txt', '')}_filter_{len(uniprot_ids)}_{name}.txt"
        
        with open(filename, 'r', encoding='utf-8') as infile, \
             open(output_path, 'w', encoding='utf-8') as outfile:
            
            # Write header
            header = infile.readline()
            outfile.write(header)
            
            # Filter lines
            for line in infile:
                line = line.strip()
                if line:
                    items = line.split('\t')
                    if items and items[0] in uniprot_ids:
                        outfile.write(f"{line}\n")
        
        self._log(f"Filtered {len(uniprot_ids)} UniProt IDs from {filename}")
        return output_path
    
    def filter_by_gene_ids(self, gene_ids: List[str], filename: str, 
                          output_path: Optional[str] = None) -> str:
        """Filter file by list of gene IDs."""
        if not gene_ids:
            raise ValueError("No gene IDs provided")
        
        if not output_path:
            name = random.choice(gene_ids)
            output_path = f"{filename.replace('.txt', '')}_filter_{len(gene_ids)}_{name}.txt"
        
        with open(filename, 'r', encoding='utf-8') as infile, \
             open(output_path, 'w', encoding='utf-8') as outfile:
            
            # Write header
            header = infile.readline()
            outfile.write(header)
            
            # Filter lines
            for line in infile:
                line = line.strip()
                if line:
                    items = line.split('\t')
                    for gene_id in gene_ids:
                        for item in items:
                            if gene_id in item.split(' '):
                                outfile.write(f"{line}\n")
                                break
        
        self._log(f"Filtered by {len(gene_ids)} gene IDs from {filename}")
        return output_path
    
    def filter_by_ipr_ids(self, ipr_ids: List[str], filename: str, 
                         output_path: Optional[str] = None) -> str:
        """Filter file by list of InterPro IDs."""
        if not ipr_ids:
            raise ValueError("No InterPro IDs provided")
        
        if not output_path:
            name = random.choice(ipr_ids)
            output_path = f"{filename.replace('.txt', '')}_filter_{len(ipr_ids)}_{name}_IPR2ACs.txt"
        
        unique_lines = set()
        
        with open(filename, 'r', encoding='utf-8') as infile, \
             open(output_path, 'w', encoding='utf-8') as outfile:
            
            # Write header
            header = infile.readline()
            outfile.write(header)
            
            # Filter lines
            for line in infile:
                line = line.strip()
                if line:
                    items = line.split('\t')
                    for item in items:
                        if item.startswith('IPR'):
                            ipr_list = item.split('; ')
                            for ipr in ipr_list:
                                if ipr in ipr_ids and line not in unique_lines:
                                    unique_lines.add(line)
                                    outfile.write(f"{line}\n")
                                    break
        
        self._log(f"Filtered by {len(ipr_ids)} InterPro IDs from {filename}")
        return output_path
    
    def filter_by_go_ids(self, go_ids: List[str], filename: str, 
                         output_path: Optional[str] = None) -> str:
        """Filter file by list of GO IDs."""
        if not go_ids:
            raise ValueError("No GO IDs provided")
        
        if not output_path:
            name = random.choice(go_ids)
            output_path = f"{filename.replace('.txt', '')}_filter_{len(go_ids)}_{name}_GO2ACs.txt"
        
        with open(filename, 'r', encoding='utf-8') as infile, \
             open(output_path, 'w', encoding='utf-8') as outfile:
            
            # Write header
            header = infile.readline()
            outfile.write(header)
            
            # Filter lines
            for line in infile:
                line = line.strip()
                if line:
                    items = line.split('\t')
                    for item in items:
                        if item.startswith('GO:'):
                            go_list = item.split('; ')
                            for go in go_list:
                                if go in go_ids:
                                    outfile.write(f"{line}\n")
                                    break
        
        self._log(f"Filtered by {len(go_ids)} GO IDs from {filename}")
        return output_path
    
    def find_words_in_file(self, words: List[str], filename: str, 
                          output_path: Optional[str] = None) -> str:
        """Find lines containing any of the specified words."""
        if not words:
            raise ValueError("No words provided")
        
        if not output_path:
            output_path = f"{filename.replace('.txt', '')}_find_{words[0]}.txt"
        
        # Expand word list with case variations
        expanded_words = []
        for word in words[:3]:  # Limit to first 3 words
            expanded_words.extend([word, word.lower(), word.upper()])
        
        unique_lines = set()
        count = 0
        
        with open(filename, 'r', encoding='utf-8') as infile, \
             open(output_path, 'w', encoding='utf-8') as outfile:
            
            # Write header
            header = infile.readline()
            outfile.write(header)
            
            # Search for words
            for line in infile:
                line = line.strip()
                if line:
                    items = line.split('\t')
                    for item in items:
                        for word in expanded_words:
                            if word in item and line not in unique_lines:
                                unique_lines.add(line)
                                outfile.write(f"{line}\n")
                                count += 1
                                break
        
        self._log(f"Found {count} lines containing words: {words}")
        return output_path
    
    def remove_duplicates(self, filename: str, output_path: Optional[str] = None) -> str:
        """Remove duplicate lines from a file."""
        if not output_path:
            output_path = f"{filename.replace('.txt', '')}_filtered.txt"
        
        unique_lines = []
        count = 0
        
        with open(filename, 'r', encoding='utf-8') as infile, \
             open(output_path, 'w', encoding='utf-8') as outfile:
            
            for line in infile:
                count += 1
                if line not in unique_lines:
                    unique_lines.append(line)
                    outfile.write(line)
        
        self._log(f"Removed {count - len(unique_lines)} duplicate lines from {filename}")
        return output_path
