"""Tests for analysis plotting functionality."""

import unittest
import tempfile
import os
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

from analysis.plotting import PlotGenerator


class TestPlotGenerator(unittest.TestCase):
    """Test cases for PlotGenerator class."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.plotter = PlotGenerator()
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """Clean up test fixtures."""
        import shutil
        shutil.rmtree(self.temp_dir)
    
    def test_create_bar_plot(self):
        """Test creating a bar plot."""
        data = {'A': 10, 'B': 20, 'C': 15}
        output_path = os.path.join(self.temp_dir, "test_bar.png")
        
        result = self.plotter.create_bar_plot(data, output_path=output_path)
        
        self.assertEqual(result, output_path)
        self.assertTrue(os.path.exists(output_path))
    
    def test_create_histogram(self):
        """Test creating a histogram."""
        data = [1.0, 2.0, 3.0, 4.0, 5.0, 2.0, 3.0, 4.0]
        output_path = os.path.join(self.temp_dir, "test_hist.png")
        
        result = self.plotter.create_histogram(data, output_path=output_path)
        
        self.assertEqual(result, output_path)
        self.assertTrue(os.path.exists(output_path))
    
    def test_create_scatter_plot(self):
        """Test creating a scatter plot."""
        x_data = [1.0, 2.0, 3.0, 4.0, 5.0]
        y_data = [2.0, 4.0, 6.0, 8.0, 10.0]
        output_path = os.path.join(self.temp_dir, "test_scatter.png")
        
        result = self.plotter.create_scatter_plot(x_data, y_data, output_path=output_path)
        
        self.assertEqual(result, output_path)
        self.assertTrue(os.path.exists(output_path))
    
    def test_create_line_plot(self):
        """Test creating a line plot."""
        x_data = [1.0, 2.0, 3.0, 4.0, 5.0]
        y_data = [2.0, 4.0, 6.0, 8.0, 10.0]
        output_path = os.path.join(self.temp_dir, "test_line.png")
        
        result = self.plotter.create_line_plot(x_data, y_data, output_path=output_path)
        
        self.assertEqual(result, output_path)
        self.assertTrue(os.path.exists(output_path))
    
    def test_save_image_png(self):
        """Test saving image in PNG format."""
        import matplotlib.pyplot as plt
        
        # Create a simple plot
        plt.figure()
        plt.plot([1, 2, 3], [1, 4, 9])
        
        output_path = os.path.join(self.temp_dir, "test_save.png")
        result = self.plotter.save_image(output_path)
        
        self.assertEqual(result, output_path)
        self.assertTrue(os.path.exists(output_path))
        plt.close()
    
    def test_save_image_pdf(self):
        """Test saving image in PDF format."""
        import matplotlib.pyplot as plt
        
        # Create a simple plot
        plt.figure()
        plt.plot([1, 2, 3], [1, 4, 9])
        
        output_path = os.path.join(self.temp_dir, "test_save.pdf")
        result = self.plotter.save_image(output_path)
        
        self.assertEqual(result, output_path)
        self.assertTrue(os.path.exists(output_path))
        plt.close()


if __name__ == '__main__':
    unittest.main()
