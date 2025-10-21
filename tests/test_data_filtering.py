"""Tests for data filtering functionality."""

import unittest
import tempfile
import os

from data.filtering import DataFilter


class TestDataFilter(unittest.TestCase):
    """Test cases for DataFilter class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.filter = DataFilter()
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
    
    def test_filter_by_uniprot_ids(self):
        """Test filtering by UniProt IDs."""
        # Create test file
        content = "ID\tName\tValue\nA1\tProtein1\t100\nA2\tProtein2\t200\nA3\tProtein3\t300\n"
        test_file = self.create_test_file("test.txt", content)
        
        # Filter by UniProt IDs
        uniprot_ids = ["A1", "A3"]
        output = self.filter.filter_by_uniprot_ids(uniprot_ids, test_file)
        
        # Check output file exists
        self.assertTrue(os.path.exists(output))
        
        # Check content
        with open(output, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 3)  # Header + 2 data lines
            self.assertTrue("A1\tProtein1\t100" in lines[1])
            self.assertTrue("A3\tProtein3\t300" in lines[2])
    
    def test_filter_by_ipr_ids(self):
        """Test filtering by InterPro IDs."""
        # Create test file
        content = "ID\tName\tIPR\nA1\tProtein1\tIPR001234; IPR005678\nA2\tProtein2\tIPR009876\nA3\tProtein3\tIPR001234\n"
        test_file = self.create_test_file("test.txt", content)
        
        # Filter by InterPro IDs
        ipr_ids = ["IPR001234"]
        output = self.filter.filter_by_ipr_ids(ipr_ids, test_file)
        
        # Check output file exists
        self.assertTrue(os.path.exists(output))
        
        # Check content
        with open(output, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 3)  # Header + 2 data lines
            self.assertTrue("A1\tProtein1\tIPR001234; IPR005678" in lines[1])
            self.assertTrue("A3\tProtein3\tIPR001234" in lines[2])
    
    def test_find_words_in_file(self):
        """Test finding words in file."""
        # Create test file
        content = "ID\tName\tDescription\nA1\tProtein1\tkinase activity\nA2\tProtein2\tbinding\nA3\tProtein3\tkinase domain\n"
        test_file = self.create_test_file("test.txt", content)
        
        # Find words
        words = ["kinase"]
        output = self.filter.find_words_in_file(words, test_file)
        
        # Check output file exists
        self.assertTrue(os.path.exists(output))
        
        # Check content
        with open(output, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 3)  # Header + 2 data lines
            self.assertTrue("A1\tProtein1\tkinase activity" in lines[1])
            self.assertTrue("A3\tProtein3\tkinase domain" in lines[2])
    
    def test_remove_duplicates(self):
        """Test removing duplicate lines."""
        # Create test file with duplicates
        content = "ID\tName\nA1\tProtein1\nA2\tProtein2\nA1\tProtein1\nA3\tProtein3\n"
        test_file = self.create_test_file("test.txt", content)
        
        # Remove duplicates
        output = self.filter.remove_duplicates(test_file)
        
        # Check output file exists
        self.assertTrue(os.path.exists(output))
        
        # Check content
        with open(output, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 4)  # Header + 3 unique data lines


if __name__ == '__main__':
    unittest.main()
