"""Tests for core functionality."""

import unittest
import tempfile
import os
from pathlib import Path

from core.config import AppConfig
from core.session import SessionManager


class TestAppConfig(unittest.TestCase):
    """Test cases for AppConfig class."""
    
    def test_default_values(self):
        """Test default configuration values."""
        config = AppConfig()
        self.assertEqual(config.app_name, "LeishDomains")
        self.assertEqual(config.default_session_prefix, "SESSION")
        self.assertIsInstance(config.cwd, str)
    
    def test_cwd_property(self):
        """Test current working directory property."""
        config = AppConfig()
        expected_cwd = os.getcwd()
        self.assertEqual(config.cwd, expected_cwd)


class TestSessionManager(unittest.TestCase):
    """Test cases for SessionManager class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()
        self.session_manager = SessionManager()
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_save_session(self):
        """Test saving a session."""
        file_paths = ["file1.txt", "file2.txt", "file3.txt"]
        output_path = self.session_manager.save(file_paths, self.temp_dir)
        
        # Check file was created
        self.assertTrue(os.path.exists(output_path))
        
        # Check content
        with open(output_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 3)
            self.assertEqual(lines[0].strip(), "file1.txt")
            self.assertEqual(lines[1].strip(), "file2.txt")
            self.assertEqual(lines[2].strip(), "file3.txt")
    
    def test_load_session(self):
        """Test loading a session."""
        # Create a session file
        session_file = os.path.join(self.temp_dir, "test_session.txt")
        with open(session_file, 'w', encoding='utf-8') as f:
            f.write("file1.txt\nfile2.txt\nfile3.txt\n")
        
        # Load session
        loaded_files = self.session_manager.load(session_file)
        
        # Check loaded files
        self.assertEqual(len(loaded_files), 3)
        self.assertEqual(loaded_files[0], "file1.txt")
        self.assertEqual(loaded_files[1], "file2.txt")
        self.assertEqual(loaded_files[2], "file3.txt")
    
    def test_save_empty_session(self):
        """Test saving empty session raises error."""
        with self.assertRaises(ValueError):
            self.session_manager.save([], self.temp_dir)
    
    def test_load_nonexistent_session(self):
        """Test loading nonexistent session raises error."""
        nonexistent_file = os.path.join(self.temp_dir, "nonexistent.txt")
        with self.assertRaises(FileNotFoundError):
            self.session_manager.load(nonexistent_file)


if __name__ == '__main__':
    unittest.main()
