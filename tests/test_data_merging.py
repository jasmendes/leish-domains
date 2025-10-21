"""Tests for data merging functionality."""

import unittest
import tempfile
import os
from pathlib import Path

from data.merging import DataMerger


class TestDataMerger(unittest.TestCase):
    """Test cases for DataMerger class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.merger = DataMerger()
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def create_test_file(self, filename: str, content: str) -> str:
        """Create a test file with given content."""
        filepath = os.path.join(self.temp_dir, filename)
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        return filepath
    
    def test_merge_by_accession(self):
        """Test merging files by accession number."""
        # Create test files
        file1_content = "ID\tName\tValue\nA1\tProtein1\t100\nA2\tProtein2\t200\n"
        file2_content = "ID\tDescription\tScore\nA1\tDesc1\t0.8\nA3\tDesc3\t0.9\n"
        
        file1 = self.create_test_file("test1.txt", file1_content)
        file2 = self.create_test_file("test2.txt", file2_content)
        
        # Merge files
        output = self.merger.merge_by_accession(file1, file2)
        
        # Check output file exists
        self.assertTrue(os.path.exists(output))
        
        # Check content
        with open(output, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 3)  # Header + 2 data lines
            self.assertTrue("A1\tProtein1\t100\tA1\tDesc1\t0.8" in lines[1])
    
    def test_match_by_accession(self):
        """Test matching files by accession number."""
        # Create test files
        file1_content = "ID\tName\tValue\nA1\tProtein1\t100\nA2\tProtein2\t200\n"
        file2_content = "ID\tDescription\tScore\nA1\tDesc1\t0.8\nA3\tDesc3\t0.9\n"
        
        file1 = self.create_test_file("test1.txt", file1_content)
        file2 = self.create_test_file("test2.txt", file2_content)
        
        # Match files
        output = self.merger.match_by_accession(file1, file2)
        
        # Check output file exists
        self.assertTrue(os.path.exists(output))
        
        # Check content (should only have A1)
        with open(output, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 2)  # Header + 1 data line
            self.assertTrue("A1\tProtein1\t100" in lines[1])
    
    def test_join_lines(self):
        """Test joining lines from two files."""
        # Create test files
        file1_content = "ID\tName\nA1\tProtein1\nA2\tProtein2\n"
        file2_content = "ID\tDescription\nA3\tDesc3\nA4\tDesc4\n"
        
        file1 = self.create_test_file("test1.txt", file1_content)
        file2 = self.create_test_file("test2.txt", file2_content)
        
        # Join files
        output = self.merger.join_lines(file1, file2)
        
        # Check output file exists
        self.assertTrue(os.path.exists(output))
        
        # Check content
        with open(output, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 5)  # Header + 4 data lines


if __name__ == '__main__':
    unittest.main()
