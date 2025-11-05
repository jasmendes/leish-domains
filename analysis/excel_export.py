"""Excel export functionality."""

import os
import xlwt
from typing import List, Optional


class ExcelExporter:
    """Handles Excel export operations."""
    
    def __init__(self):
        self.logger = None  # Will be set by caller if needed
    
    def _log(self, message: str):
        """Log message if logger is available."""
        if self.logger:
            self.logger.info(message)
        else:
            print(message)
    
    def _is_number(self, s: str) -> bool:
        """Check if a string represents a number."""
        try:
            float(s)
            return True
        except ValueError:
            return False
    
    def convert_txt_to_excel(self, textfile: str, output_path: Optional[str] = None) -> str:
        """Convert a text file to Excel format."""
        if not output_path:
            output_path = textfile.replace('.txt', '.xls')
        
        # Read text file
        with open(textfile, 'r', encoding='utf-8') as f:
            row_list = []
            for row in f:
                row_list.append(row.split('\t'))
        
        # Create Excel workbook
        workbook = xlwt.Workbook()
        worksheet = workbook.add_sheet('Sheet1')
        
        # Write data to Excel
        for row_idx, row in enumerate(row_list):
            for col_idx, cell in enumerate(row):
                value = cell.strip()
                if self._is_number(value):
                    try:
                        worksheet.write(row_idx, col_idx, float(value))
                    except ValueError:
                        worksheet.write(row_idx, col_idx, value)
                else:
                    worksheet.write(row_idx, col_idx, value)
        
        # Save workbook
        workbook.save(output_path)
        self._log(f"Converted {textfile} to {output_path}")
        
        return output_path
    
    def convert_multiple_files_to_excel(self, file_paths: List[str], 
                                      output_path: str) -> str:
        """Convert multiple text files to a single Excel workbook."""
        workbook = xlwt.Workbook()
        
        for sheet_idx, file_path in enumerate(file_paths):
            sheet_name = f"Sheet{sheet_idx}"
            worksheet = workbook.add_sheet(sheet_name)
            
            # Read and write file data
            with open(file_path, 'r', encoding='utf-8') as f:
                row_list = []
                for row in f:
                    row_list.append(row.split('\t'))
                
                for row_idx, row in enumerate(row_list):
                    for col_idx, cell in enumerate(row):
                        value = cell.strip()
                        if self._is_number(value):
                            try:
                                worksheet.write(row_idx, col_idx, float(value))
                            except ValueError:
                                worksheet.write(row_idx, col_idx, value)
                        else:
                            worksheet.write(row_idx, col_idx, value)
            
            # Add file info
            worksheet.write(0, len(row_list[0]) if row_list else 0, "File_info")
            worksheet.write(1, len(row_list[0]) if row_list else 0, file_path)
        
        # Save workbook
        workbook.save(output_path)
        self._log(f"Converted {len(file_paths)} files to {output_path}")
        
        return output_path
    
    def create_excel_worksheet(self, directory_path: str, delimiter: str = '|') -> None:
        """Create Excel worksheets for all text files in a directory."""
        from os import listdir
        from os.path import isfile, join
        
        # Get all text files
        textfiles = [join(directory_path, f) for f in listdir(directory_path) 
                    if isfile(join(directory_path, f)) and f.endswith('.txt')]
        
        for textfile in textfiles:
            try:
                self.convert_txt_to_excel(textfile)
            except Exception as e:
                self._log(f"Error converting {textfile}: {e}")
        
        self._log(f"Converted {len(textfiles)} files in {directory_path}")
