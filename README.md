# LeishDomains

A comprehensive Python application for analyzing protein domains in Leishmania species. The application provides both a graphical user interface (GUI) and a command-line interface (CLI) for flexible data analysis workflows.

## Features

### Core Functionality
- **Session Management**: Create, load, and save analysis sessions
- **Data Merging**: Merge files by accession numbers, InterPro IDs, or other criteria
- **Data Filtering**: Filter datasets by UniProt IDs, gene IDs, InterPro IDs, GO terms, or keywords
- **Statistical Analysis**: Calculate statistics and generate visualizations
- **Excel Export**: Convert text files to Excel format with multiple sheets

### User Interfaces
- **GUI Mode**: User-friendly Tkinter interface for interactive analysis
- **CLI Mode**: Command-line interface for batch processing and automation
- **Dual Mode**: Seamless switching between GUI and CLI

### Data Processing
- File operations (slicing, cutting, splitting)
- Data merging and joining
- Duplicate removal
- Text search and filtering
- Statistical analysis and plotting

## Installation

1. Clone this repository:
```bash
git clone <repository-url>
cd leish-domains
```

2. Install the required dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### GUI Mode (Default)
```bash
python main.py
# or
python main.py --gui
```

### CLI Mode
```bash
python main.py --cli <command> [options]
```

### CLI Commands

#### Merge Files
```bash
# Merge two files by accession number
python main.py merge file1.txt file2.txt -o merged.txt

# Match only (keep only entries present in both files)
python main.py merge file1.txt file2.txt --match-only -o matched.txt
```

#### Filter Data
```bash
# Filter by UniProt IDs
python main.py filter data.txt --uniprot-ids A1B2C3,D4E5F6 -o filtered.txt

# Filter by InterPro IDs
python main.py filter data.txt --ipr-ids IPR001234,IPR005678 -o filtered.txt

# Filter by GO terms
python main.py filter data.txt --go-ids GO:0003674,GO:0008150 -o filtered.txt

# Search for keywords
python main.py filter data.txt --words kinase,phosphatase -o filtered.txt
```

#### Analyze Data
```bash
# Analyze a specific column
python main.py analyze data.txt --column 2

# Generate histogram
python main.py analyze data.txt --column 2 --plot histogram.png --plot-type hist
```

#### Session Management
```bash
# Save session
python main.py session save --files file1.txt,file2.txt,file3.txt -o session.txt

# Load session
python main.py session load session.txt
```

#### Excel Export
```bash
# Convert single file to Excel
python main.py excel data.txt -o output.xls

# Convert multiple files to Excel workbook
python main.py excel file1.txt file2.txt file3.txt -o workbook.xls
```

## Project Structure

```
leish-domains/
├── core/                   # Core utilities
│   ├── __init__.py
│   ├── config.py          # Configuration management
│   ├── logger.py          # Logging utilities
│   └── session.py         # Session management
├── data/                   # Data processing modules
│   ├── __init__.py
│   ├── file_ops.py        # File operations
│   ├── merging.py         # Data merging
│   └── filtering.py        # Data filtering
├── analysis/              # Analysis modules
│   ├── __init__.py
│   ├── statistics.py      # Statistical analysis
│   ├── plotting.py        # Plotting utilities
│   └── excel_export.py    # Excel export
├── cli/                   # Command-line interface
│   ├── __init__.py
│   └── main.py           # CLI implementation
├── tests/                 # Test suite
│   ├── __init__.py
│   ├── test_core.py
│   ├── test_data_merging.py
│   └── test_data_filtering.py
├── main.py               # Main application entry point
├── requirements.txt      # Dependencies
└── README.md            # This file

leish-domains/
├── .github/workflows/     # CI/CD pipeline
├── core/                  # Core utilities + config + validation
├── data/                  # Data processing modules
├── analysis/              # Analysis modules
├── cli/                   # Command-line interface
├── tests/                 # Test suite
├── examples/              # Usage examples and sample data
├── docs/                  # Sphinx documentation
├── config/                # Configuration files
├── main.py               # Dual-mode entry point
├── setup.py              # Package setup
├── pyproject.toml        # Modern Python packaging
├── Makefile              # Development tasks
├── requirements.txt      # Dependencies
└── README.md             # Comprehensive documentation
```

## Dependencies

- Python 3.8+
- Tkinter (usually comes with Python)
- Biopython >= 1.79
- NumPy >= 1.24.0
- Matplotlib >= 3.7.0
- SciPy >= 1.10.0
- xlwt (for Excel export)

## Testing

Run the test suite:
```bash
python -m pytest tests/
```

Or run specific test modules:
```bash
python -m pytest tests/test_core.py
python -m pytest tests/test_data_merging.py
python -m pytest tests/test_data_filtering.py
```

## Development

### Adding New Features
1. Create new modules in appropriate directories (`data/`, `analysis/`, etc.)
2. Add corresponding tests in `tests/`
3. Update CLI interface in `cli/main.py` if needed
4. Update documentation

### Code Style
- Follow PEP 8 guidelines
- Use type hints where appropriate
- Add docstrings to all public methods
- Write tests for new functionality

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## Support

For issues and questions, please create an issue in the repository or contact the development team.