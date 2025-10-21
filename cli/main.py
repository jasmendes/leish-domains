"""Command-line interface implementation."""

import argparse
import sys
import os
from typing import List, Optional
from pathlib import Path

from core.logger import get_logger
from core.session import SessionManager
from data import FileOperations, DataMerger, DataFilter
from analysis import DataAnalyzer, PlotGenerator, ExcelExporter


class CLIApp:
    """Command-line interface for LeishDomains."""
    
    def __init__(self):
        self.logger = get_logger("leishdomains.cli")
        self.session_manager = SessionManager()
        self.file_ops = FileOperations()
        self.data_merger = DataMerger()
        self.data_filter = DataFilter()
        self.analyzer = DataAnalyzer()
        self.plotter = PlotGenerator()
        self.excel_exporter = ExcelExporter()
        
        # Set loggers for sub-modules
        self.data_merger.logger = self.logger
        self.data_filter.logger = self.logger
        self.analyzer.logger = self.logger
        self.plotter.logger = self.logger
        self.excel_exporter.logger = self.logger
    
    def create_parser(self) -> argparse.ArgumentParser:
        """Create command-line argument parser."""
        parser = argparse.ArgumentParser(
            description="LeishDomains - Protein domain analysis tool",
            formatter_class=argparse.RawDescriptionHelpFormatter,
            epilog="""
Examples:
  %(prog)s merge file1.txt file2.txt -o merged.txt
  %(prog)s filter --uniprot-ids A1B2C3,D4E5F6 data.txt
  %(prog)s analyze data.txt --column 2 --plot histogram.png
  %(prog)s session save --files file1.txt,file2.txt
  %(prog)s session load session.txt
            """
        )
        
        subparsers = parser.add_subparsers(dest='command', help='Available commands')
        
        # Merge command
        merge_parser = subparsers.add_parser('merge', help='Merge two files')
        merge_parser.add_argument('file1', help='First file to merge')
        merge_parser.add_argument('file2', help='Second file to merge')
        merge_parser.add_argument('-o', '--output', help='Output file path')
        merge_parser.add_argument('--match-only', action='store_true', 
                                 help='Only keep matching entries')
        
        # Filter command
        filter_parser = subparsers.add_parser('filter', help='Filter data by criteria')
        filter_parser.add_argument('input_file', help='Input file to filter')
        filter_parser.add_argument('-o', '--output', help='Output file path')
        filter_parser.add_argument('--uniprot-ids', help='Comma-separated UniProt IDs')
        filter_parser.add_argument('--gene-ids', help='Comma-separated gene IDs')
        filter_parser.add_argument('--ipr-ids', help='Comma-separated InterPro IDs')
        filter_parser.add_argument('--go-ids', help='Comma-separated GO IDs')
        filter_parser.add_argument('--words', help='Comma-separated words to search for')
        
        # Analyze command
        analyze_parser = subparsers.add_parser('analyze', help='Analyze data')
        analyze_parser.add_argument('input_file', help='Input file to analyze')
        analyze_parser.add_argument('--column', type=int, help='Column index to analyze')
        analyze_parser.add_argument('--plot', help='Generate plot and save to file')
        analyze_parser.add_argument('--plot-type', choices=['bar', 'hist', 'scatter', 'line'],
                                   default='hist', help='Type of plot to generate')
        
        # Session commands
        session_parser = subparsers.add_parser('session', help='Session management')
        session_subparsers = session_parser.add_subparsers(dest='session_action')
        
        # Save session
        save_parser = session_subparsers.add_parser('save', help='Save session')
        save_parser.add_argument('--files', required=True, 
                               help='Comma-separated list of files')
        save_parser.add_argument('-o', '--output', help='Output session file')
        
        # Load session
        load_parser = session_subparsers.add_parser('load', help='Load session')
        load_parser.add_argument('session_file', help='Session file to load')
        
        # Excel export
        excel_parser = subparsers.add_parser('excel', help='Export to Excel')
        excel_parser.add_argument('input_files', nargs='+', help='Input files to convert')
        excel_parser.add_argument('-o', '--output', help='Output Excel file')
        
        return parser
    
    def run_merge(self, args) -> int:
        """Run merge command."""
        try:
            if args.match_only:
                output = self.data_merger.match_by_accession(args.file1, args.file2, args.output)
            else:
                output = self.data_merger.merge_by_accession(args.file1, args.file2, args.output)
            
            self.logger.info(f"Merge completed: {output}")
            return 0
        except Exception as e:
            self.logger.error(f"Merge failed: {e}")
            return 1
    
    def run_filter(self, args) -> int:
        """Run filter command."""
        try:
            output = None
            
            if args.uniprot_ids:
                ids = [id.strip() for id in args.uniprot_ids.split(',')]
                output = self.data_filter.filter_by_uniprot_ids(ids, args.input_file, args.output)
            elif args.gene_ids:
                ids = [id.strip() for id in args.gene_ids.split(',')]
                output = self.data_filter.filter_by_gene_ids(ids, args.input_file, args.output)
            elif args.ipr_ids:
                ids = [id.strip() for id in args.ipr_ids.split(',')]
                output = self.data_filter.filter_by_ipr_ids(ids, args.input_file, args.output)
            elif args.go_ids:
                ids = [id.strip() for id in args.go_ids.split(',')]
                output = self.data_filter.filter_by_go_ids(ids, args.input_file, args.output)
            elif args.words:
                words = [word.strip() for word in args.words.split(',')]
                output = self.data_filter.find_words_in_file(words, args.input_file, args.output)
            else:
                self.logger.error("No filter criteria specified")
                return 1
            
            self.logger.info(f"Filter completed: {output}")
            return 0
        except Exception as e:
            self.logger.error(f"Filter failed: {e}")
            return 1
    
    def run_analyze(self, args) -> int:
        """Run analyze command."""
        try:
            if args.column is not None:
                stats = self.analyzer.analyze_column(args.input_file, args.column)
                self.logger.info(f"Column {args.column} statistics: {stats}")
            
            if args.plot:
                # Generate plot based on type
                if args.plot_type == 'hist':
                    # For histogram, we need numeric data
                    with open(args.input_file, 'r') as f:
                        f.readline()  # Skip header
                        data = []
                        for line in f:
                            items = line.strip().split('\t')
                            if len(items) > args.column and self.analyzer.is_number(items[args.column]):
                                data.append(float(items[args.column]))
                    
                    if data:
                        self.plotter.create_histogram(data, output_path=args.plot)
                        self.logger.info(f"Histogram saved: {args.plot}")
                    else:
                        self.logger.error("No numeric data found for histogram")
                        return 1
                else:
                    self.logger.error(f"Plot type {args.plot_type} not implemented yet")
                    return 1
            
            return 0
        except Exception as e:
            self.logger.error(f"Analysis failed: {e}")
            return 1
    
    def run_session(self, args) -> int:
        """Run session command."""
        try:
            if args.session_action == 'save':
                files = [f.strip() for f in args.files.split(',')]
                output = self.session_manager.save(files)
                self.logger.info(f"Session saved: {output}")
                return 0
            elif args.session_action == 'load':
                files = self.session_manager.load(args.session_file)
                self.logger.info(f"Session loaded: {len(files)} files")
                for file in files:
                    print(file)
                return 0
            else:
                self.logger.error("No session action specified")
                return 1
        except Exception as e:
            self.logger.error(f"Session operation failed: {e}")
            return 1
    
    def run_excel(self, args) -> int:
        """Run Excel export command."""
        try:
            if len(args.input_files) == 1:
                output = self.excel_exporter.convert_txt_to_excel(args.input_files[0], args.output)
            else:
                output = self.excel_exporter.convert_multiple_files_to_excel(args.input_files, args.output)
            
            self.logger.info(f"Excel export completed: {output}")
            return 0
        except Exception as e:
            self.logger.error(f"Excel export failed: {e}")
            return 1
    
    def run(self, args: Optional[List[str]] = None) -> int:
        """Run the CLI application."""
        parser = self.create_parser()
        parsed_args = parser.parse_args(args)
        
        if not parsed_args.command:
            parser.print_help()
            return 1
        
        # Route to appropriate handler
        if parsed_args.command == 'merge':
            return self.run_merge(parsed_args)
        elif parsed_args.command == 'filter':
            return self.run_filter(parsed_args)
        elif parsed_args.command == 'analyze':
            return self.run_analyze(parsed_args)
        elif parsed_args.command == 'session':
            return self.run_session(parsed_args)
        elif parsed_args.command == 'excel':
            return self.run_excel(parsed_args)
        else:
            self.logger.error(f"Unknown command: {parsed_args.command}")
            return 1


def main():
    """Main entry point for CLI."""
    app = CLIApp()
    sys.exit(app.run())


if __name__ == "__main__":
    main()
