"""File operations and utilities."""

import os
import itertools
from typing import List, Iterator, Optional
from pathlib import Path


class FileOperations:
    """Handles basic file operations and utilities."""
    
    @staticmethod
    def slice_file(filename: str, start: int, end: int) -> Iterator[str]:
        """Slice a file to get lines from start to end (exclusive)."""
        with open(filename, 'r', encoding='utf-8') as f:
            return itertools.islice(f, start, end)
    
    @staticmethod
    def cut_file(filename: str, start: int, end: int, output_path: str = "output/new_file.txt") -> str:
        """Cut lines from a file and save to output."""
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        with open(output_path, 'w', encoding='utf-8') as out:
            for line in FileOperations.slice_file(filename, start, end):
                out.write(line)
        
        return output_path
    
    @staticmethod
    def split_file_by_delimiter(filename: str, delimiter: str = "kinase\n") -> tuple[str, str]:
        """Split a file by delimiter and return two parts."""
        with open(filename, 'r', encoding='utf-8') as f:
            content = f.read()
        
        part1, sentinel, part2 = content.partition(delimiter)
        
        part1_path = filename.replace('.txt', '_p1.txt')
        part2_path = filename.replace('.txt', '_p2.txt')
        
        with open(part1_path, 'w', encoding='utf-8') as f:
            f.write(part1)
        with open(part2_path, 'w', encoding='utf-8') as f:
            f.write(part2)
        
        return part1_path, part2_path
    
    @staticmethod
    def is_number(s: str) -> bool:
        """Check if a string represents a number."""
        try:
            float(s)
            return True
        except ValueError:
            return False
    
    @staticmethod
    def remove_exponent(value) -> str:
        """Remove scientific notation from decimal values."""
        import decimal
        
        decimal_places = 8
        max_digits = 16
        
        if isinstance(value, decimal.Decimal):
            context = decimal.getcontext().copy()
            context.prec = max_digits
            return "{0:f}".format(value.quantize(decimal.Decimal(".1") ** decimal_places, context=context))
        else:
            return "%.*f" % (decimal_places, value)
