"""Plotting and visualization utilities."""

import matplotlib.pyplot as plt
import numpy as np
from typing import List, Dict, Optional, Tuple
from pathlib import Path


class PlotGenerator:
    """Handles plotting and visualization of data."""
    
    def __init__(self):
        self.logger = None  # Will be set by caller if needed
    
    def _log(self, message: str):
        """Log message if logger is available."""
        if self.logger:
            self.logger.info(message)
        else:
            print(message)
    
    def save_image(self, output_path: str, file_types: List[str] = None) -> str:
        """Save current plot to file with specified format."""
        if file_types is None:
            file_types = ['PNG', 'JPEG', 'PDF']
        
        # Determine format from file extension
        ext = Path(output_path).suffix.lower()
        if ext == '.png':
            plt.savefig(output_path, format='png', dpi=300, bbox_inches='tight')
        elif ext == '.jpg' or ext == '.jpeg':
            plt.savefig(output_path, format='jpeg', dpi=300, bbox_inches='tight')
        elif ext == '.pdf':
            plt.savefig(output_path, format='pdf', bbox_inches='tight')
        else:
            plt.savefig(output_path, dpi=300, bbox_inches='tight')
        
        self._log(f"Plot saved to: {output_path}")
        return output_path
    
    def create_bar_plot(self, data: Dict[str, float], title: str = "Bar Plot", 
                       xlabel: str = "Categories", ylabel: str = "Values",
                       output_path: Optional[str] = None) -> str:
        """Create a bar plot from data dictionary."""
        plt.figure(figsize=(10, 6))
        
        categories = list(data.keys())
        values = list(data.values())
        
        plt.bar(categories, values)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        if output_path:
            self.save_image(output_path)
        
        return output_path or "bar_plot.png"
    
    def create_histogram(self, data: List[float], title: str = "Histogram",
                        bins: int = 30, output_path: Optional[str] = None) -> str:
        """Create a histogram from data."""
        plt.figure(figsize=(10, 6))
        
        plt.hist(data, bins=bins, alpha=0.7, edgecolor='black')
        plt.title(title)
        plt.xlabel("Value")
        plt.ylabel("Frequency")
        plt.tight_layout()
        
        if output_path:
            self.save_image(output_path)
        
        return output_path or "histogram.png"
    
    def create_scatter_plot(self, x_data: List[float], y_data: List[float],
                           title: str = "Scatter Plot", xlabel: str = "X",
                           ylabel: str = "Y", output_path: Optional[str] = None) -> str:
        """Create a scatter plot from x and y data."""
        plt.figure(figsize=(10, 6))
        
        plt.scatter(x_data, y_data, alpha=0.6)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        
        if output_path:
            self.save_image(output_path)
        
        return output_path or "scatter_plot.png"
    
    def create_line_plot(self, x_data: List[float], y_data: List[float],
                        title: str = "Line Plot", xlabel: str = "X",
                        ylabel: str = "Y", output_path: Optional[str] = None) -> str:
        """Create a line plot from x and y data."""
        plt.figure(figsize=(10, 6))
        
        plt.plot(x_data, y_data, marker='o', linewidth=2, markersize=4)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        if output_path:
            self.save_image(output_path)
        
        return output_path or "line_plot.png"
