"""Tests for CLI functionality."""

import unittest
import tempfile
import os
import sys
from unittest.mock import patch, MagicMock

from cli.main import CLIApp


class TestCLIApp(unittest.TestCase):
    """Test cases for CLIApp class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.cli_app = CLIApp()
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
    
    def test_create_parser(self):
        """Test creating argument parser."""
        parser = self.cli_app.create_parser()
        self.assertIsNotNone(parser)
        
        # Test that parser has expected subcommands
        subcommands = ['merge', 'filter', 'analyze', 'session', 'excel']
        for cmd in subcommands:
            self.assertIn(cmd, parser._subparsers._name_parser_map)
    
    @patch('cli.main.CLIApp.run_merge')
    def test_run_merge_command(self, mock_run_merge):
        """Test running merge command."""
        mock_run_merge.return_value = 0
        
        # Create test files
        file1 = self.create_test_file("test1.txt", "ID\tName\nA1\tProtein1\n")
        file2 = self.create_test_file("test2.txt", "ID\tValue\nA1\t100\n")
        
        # Mock argparse to simulate command line arguments
        with patch('sys.argv', ['cli', 'merge', file1, file2]):
            result = self.cli_app.run(['merge', file1, file2])
            self.assertEqual(result, 0)
            mock_run_merge.assert_called_once()
    
    @patch('cli.main.CLIApp.run_filter')
    def test_run_filter_command(self, mock_run_filter):
        """Test running filter command."""
        mock_run_filter.return_value = 0
        
        # Create test file
        test_file = self.create_test_file("test.txt", "ID\tName\nA1\tProtein1\n")
        
        # Mock argparse to simulate command line arguments
        with patch('sys.argv', ['cli', 'filter', test_file, '--uniprot-ids', 'A1']):
            result = self.cli_app.run(['filter', test_file, '--uniprot-ids', 'A1'])
            self.assertEqual(result, 0)
            mock_run_filter.assert_called_once()
    
    @patch('cli.main.CLIApp.run_analyze')
    def test_run_analyze_command(self, mock_run_analyze):
        """Test running analyze command."""
        mock_run_analyze.return_value = 0
        
        # Create test file
        test_file = self.create_test_file("test.txt", "ID\tName\tValue\nA1\tProtein1\t100\n")
        
        # Mock argparse to simulate command line arguments
        with patch('sys.argv', ['cli', 'analyze', test_file, '--column', '2']):
            result = self.cli_app.run(['analyze', test_file, '--column', '2'])
            self.assertEqual(result, 0)
            mock_run_analyze.assert_called_once()
    
    @patch('cli.main.CLIApp.run_session')
    def test_run_session_command(self, mock_run_session):
        """Test running session command."""
        mock_run_session.return_value = 0
        
        # Mock argparse to simulate command line arguments
        with patch('sys.argv', ['cli', 'session', 'save', '--files', 'file1.txt,file2.txt']):
            result = self.cli_app.run(['session', 'save', '--files', 'file1.txt,file2.txt'])
            self.assertEqual(result, 0)
            mock_run_session.assert_called_once()
    
    @patch('cli.main.CLIApp.run_excel')
    def test_run_excel_command(self, mock_run_excel):
        """Test running excel command."""
        mock_run_excel.return_value = 0
        
        # Create test file
        test_file = self.create_test_file("test.txt", "ID\tName\nA1\tProtein1\n")
        
        # Mock argparse to simulate command line arguments
        with patch('sys.argv', ['cli', 'excel', test_file]):
            result = self.cli_app.run(['excel', test_file])
            self.assertEqual(result, 0)
            mock_run_excel.assert_called_once()
    
    def test_run_unknown_command(self):
        """Test running unknown command."""
        result = self.cli_app.run(['unknown_command'])
        self.assertEqual(result, 1)
    
    def test_run_merge_success(self):
        """Test successful merge operation."""
        # Create test files
        file1 = self.create_test_file("test1.txt", "ID\tName\nA1\tProtein1\n")
        file2 = self.create_test_file("test2.txt", "ID\tValue\nA1\t100\n")
        
        # Mock the data merger
        with patch.object(self.cli_app.data_merger, 'merge_by_accession') as mock_merge:
            mock_merge.return_value = "merged.txt"
            
            # Create mock args
            args = MagicMock()
            args.file1 = file1
            args.file2 = file2
            args.output = None
            args.match_only = False
            
            result = self.cli_app.run_merge(args)
            self.assertEqual(result, 0)
            mock_merge.assert_called_once_with(file1, file2, None)
    
    def test_run_filter_success(self):
        """Test successful filter operation."""
        # Create test file
        test_file = self.create_test_file("test.txt", "ID\tName\nA1\tProtein1\n")
        
        # Mock the data filter
        with patch.object(self.cli_app.data_filter, 'filter_by_uniprot_ids') as mock_filter:
            mock_filter.return_value = "filtered.txt"
            
            # Create mock args
            args = MagicMock()
            args.input_file = test_file
            args.output = None
            args.uniprot_ids = "A1,A2"
            args.gene_ids = None
            args.ipr_ids = None
            args.go_ids = None
            args.words = None
            
            result = self.cli_app.run_filter(args)
            self.assertEqual(result, 0)
            mock_filter.assert_called_once_with(['A1', 'A2'], test_file, None)
    
    def test_run_filter_no_criteria(self):
        """Test filter operation with no criteria."""
        # Create test file
        test_file = self.create_test_file("test.txt", "ID\tName\nA1\tProtein1\n")
        
        # Create mock args with no filter criteria
        args = MagicMock()
        args.input_file = test_file
        args.output = None
        args.uniprot_ids = None
        args.gene_ids = None
        args.ipr_ids = None
        args.go_ids = None
        args.words = None
        
        result = self.cli_app.run_filter(args)
        self.assertEqual(result, 1)  # Should return error code
    
    def test_run_session_save(self):
        """Test session save operation."""
        # Mock the session manager
        with patch.object(self.cli_app.session_manager, 'save') as mock_save:
            mock_save.return_value = "session.txt"
            
            # Create mock args
            args = MagicMock()
            args.session_action = 'save'
            args.files = "file1.txt,file2.txt"
            args.output = None
            
            result = self.cli_app.run_session(args)
            self.assertEqual(result, 0)
            mock_save.assert_called_once_with(['file1.txt', 'file2.txt'])
    
    def test_run_session_load(self):
        """Test session load operation."""
        # Create test session file
        session_file = self.create_test_file("session.txt", "file1.txt\nfile2.txt\n")
        
        # Mock the session manager
        with patch.object(self.cli_app.session_manager, 'load') as mock_load:
            mock_load.return_value = ['file1.txt', 'file2.txt']
            
            # Create mock args
            args = MagicMock()
            args.session_action = 'load'
            args.session_file = session_file
            
            result = self.cli_app.run_session(args)
            self.assertEqual(result, 0)
            mock_load.assert_called_once_with(session_file)


if __name__ == '__main__':
    unittest.main()
