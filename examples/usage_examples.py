#!/usr/bin/env python3
"""
Usage examples for LeishDomains CLI and programmatic API.
"""

import os
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from core.session import SessionManager
from data.merging import DataMerger
from data.filtering import DataFilter
from analysis.statistics import DataAnalyzer
from analysis.plotting import PlotGenerator
from analysis.excel_export import ExcelExporter


def example_merge_files():
    """Example: Merge two data files."""
    print("=== Merging Files Example ===")
    
    # Create sample data
    data1 = """ID\tName\tValue
A1\tProtein1\t100
A2\tProtein2\t200
A3\tProtein3\t300"""
    
    data2 = """ID\tDescription\tScore
A1\tKinase domain\t0.8
A2\tBinding domain\t0.9
A4\tRegulatory domain\t0.7"""
    
    # Write sample files
    with open("temp_file1.txt", "w") as f:
        f.write(data1)
    with open("temp_file2.txt", "w") as f:
        f.write(data2)
    
    # Merge files
    merger = DataMerger()
    output = merger.merge_by_accession("temp_file1.txt", "temp_file2.txt", "merged_output.txt")
    print(f"Merged files saved to: {output}")
    
    # Cleanup
    os.remove("temp_file1.txt")
    os.remove("temp_file2.txt")
    os.remove("merged_output.txt")


def example_filter_data():
    """Example: Filter data by criteria."""
    print("\n=== Filtering Data Example ===")
    
    # Create sample data
    data = """ID\tName\tIPR\tGO
A1\tProtein1\tIPR001234; IPR005678\tGO:0003674; GO:0008150
A2\tProtein2\tIPR009876\tGO:0003674
A3\tProtein3\tIPR001234\tGO:0008150"""
    
    with open("temp_data.txt", "w") as f:
        f.write(data)
    
    # Filter by InterPro IDs
    filter_tool = DataFilter()
    output = filter_tool.filter_by_ipr_ids(["IPR001234"], "temp_data.txt", "filtered_output.txt")
    print(f"Filtered data saved to: {output}")
    
    # Cleanup
    os.remove("temp_data.txt")
    os.remove("filtered_output.txt")


def example_analyze_data():
    """Example: Analyze data and generate plots."""
    print("\n=== Data Analysis Example ===")
    
    # Create sample data with numeric values
    data = """ID\tName\tValue
A1\tProtein1\t100.5
A2\tProtein2\t75.2
A3\tProtein3\t200.1
A4\tProtein4\t150.8
A5\tProtein5\t89.3"""
    
    with open("temp_analysis.txt", "w") as f:
        f.write(data)
    
    # Analyze data
    analyzer = DataAnalyzer()
    stats = analyzer.analyze_column("temp_analysis.txt", 2)  # Column 2 (Value)
    print(f"Statistics: {stats}")
    
    # Generate plot
    plotter = PlotGenerator()
    plotter.create_histogram([100.5, 75.2, 200.1, 150.8, 89.3], 
                            title="Protein Value Distribution",
                            output_path="protein_histogram.png")
    print("Histogram saved to: protein_histogram.png")
    
    # Cleanup
    os.remove("temp_analysis.txt")
    os.remove("protein_histogram.png")


def example_session_management():
    """Example: Session management."""
    print("\n=== Session Management Example ===")
    
    # Create sample files
    files = ["file1.txt", "file2.txt", "file3.txt"]
    for file in files:
        with open(file, "w") as f:
            f.write(f"Sample data for {file}")
    
    # Save session
    session_manager = SessionManager()
    session_path = session_manager.save(files)
    print(f"Session saved to: {session_path}")
    
    # Load session
    loaded_files = session_manager.load(session_path)
    print(f"Loaded {len(loaded_files)} files from session")
    
    # Cleanup
    for file in files:
        os.remove(file)
    os.remove(session_path)


def example_excel_export():
    """Example: Export to Excel."""
    print("\n=== Excel Export Example ===")
    
    # Create sample data
    data = """ID\tName\tValue
A1\tProtein1\t100
A2\tProtein2\t200
A3\tProtein3\t300"""
    
    with open("temp_excel.txt", "w") as f:
        f.write(data)
    
    # Export to Excel
    excel_exporter = ExcelExporter()
    output = excel_exporter.convert_txt_to_excel("temp_excel.txt", "sample_output.xls")
    print(f"Excel file saved to: {output}")
    
    # Cleanup
    os.remove("temp_excel.txt")
    os.remove("sample_output.xls")


def main():
    """Run all examples."""
    print("LeishDomains Usage Examples")
    print("=" * 40)
    
    try:
        example_merge_files()
        example_filter_data()
        example_analyze_data()
        example_session_management()
        example_excel_export()
        
        print("\n" + "=" * 40)
        print("All examples completed successfully!")
        
    except Exception as e:
        print(f"Error running examples: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
