"""Tests for analysis statistics functionality."""

import unittest
import tempfile
import os
import decimal

from analysis.statistics import DataAnalyzer


class TestDataAnalyzer(unittest.TestCase):
    """Test cases for DataAnalyzer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = DataAnalyzer()
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
    
    def test_count_frequencies(self):
        """Test counting frequencies of items."""
        data = ['a', 'b', 'a', 'c', 'b', 'a']
        result = self.analyzer.count_frequencies(data)
        expected = {'a': 3, 'b': 2, 'c': 1}
        self.assertEqual(result, expected)
    
    def test_find_most_common(self):
        """Test finding most common items."""
        data = ['a', 'b', 'a', 'c', 'b', 'a', 'd']
        result = self.analyzer.find_most_common(data, n=3)
        expected = [('a', 3), ('b', 2), ('c', 1)]
        self.assertEqual(result, expected)
    
    def test_calculate_statistics(self):
        """Test calculating basic statistics."""
        values = [1.0, 2.0, 3.0, 4.0, 5.0]
        result = self.analyzer.calculate_statistics(values)
        
        self.assertEqual(result['count'], 5)
        self.assertEqual(result['mean'], 3.0)
        self.assertEqual(result['median'], 3.0)
        self.assertEqual(result['min'], 1.0)
        self.assertEqual(result['max'], 5.0)
        self.assertEqual(result['sum'], 15.0)
    
    def test_calculate_statistics_empty(self):
        """Test calculating statistics with empty data."""
        result = self.analyzer.calculate_statistics([])
        self.assertEqual(result, {})
    
    def test_remove_exponent_decimal(self):
        """Test removing exponent from decimal values."""
        value = decimal.Decimal('1.23E+2')
        result = self.analyzer.remove_exponent(value)
        self.assertEqual(result, '123.00000000')
    
    def test_remove_exponent_float(self):
        """Test removing exponent from float values."""
        value = 1.23E+2
        result = self.analyzer.remove_exponent(value)
        self.assertEqual(result, '123.00000000')
    
    def test_is_number(self):
        """Test checking if string is a number."""
        self.assertTrue(self.analyzer.is_number('123'))
        self.assertTrue(self.analyzer.is_number('123.45'))
        self.assertTrue(self.analyzer.is_number('-123.45'))
        self.assertFalse(self.analyzer.is_number('abc'))
        self.assertFalse(self.analyzer.is_number(''))
    
    def test_analyze_column(self):
        """Test analyzing a specific column in a file."""
        content = "ID\tName\tValue\nA1\tProtein1\t100\nA2\tProtein2\t200\nA3\tProtein3\tabc\n"
        test_file = self.create_test_file("test.txt", content)
        
        result = self.analyzer.analyze_column(test_file, 2)  # Column 2 (Value)
        
        self.assertEqual(result['count'], 2)
        self.assertEqual(result['mean'], 150.0)
        self.assertEqual(result['median'], 150.0)
        self.assertEqual(result['min'], 100.0)
        self.assertEqual(result['max'], 200.0)
        self.assertEqual(result['non_numeric_count'], 1)
    
    def test_analyze_column_no_numeric(self):
        """Test analyzing column with no numeric values."""
        content = "ID\tName\tValue\nA1\tProtein1\tabc\nA2\tProtein2\tdef\n"
        test_file = self.create_test_file("test.txt", content)
        
        result = self.analyzer.analyze_column(test_file, 2)
        
        self.assertEqual(result, {'error': 'No numeric values found'})


if __name__ == '__main__':
    unittest.main()
