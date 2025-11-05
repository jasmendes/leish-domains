"""Data validation utilities."""

import os
from typing import List

from .exceptions import ValidationError, FileOperationError


class DataValidator:
    """Validates data and file operations."""
    
    @staticmethod
    def validate_file_exists(file_path: str) -> None:
        """Validate that a file exists."""
        if not os.path.isfile(file_path):
            raise FileOperationError(f"File not found: {file_path}")
    
    @staticmethod
    def validate_file_readable(file_path: str) -> None:
        """Validate that a file is readable."""
        DataValidator.validate_file_exists(file_path)
        if not os.access(file_path, os.R_OK):
            raise FileOperationError(f"File not readable: {file_path}")
    
    @staticmethod
    def validate_file_writable(file_path: str) -> None:
        """Validate that a file path is writable."""
        directory = os.path.dirname(file_path)
        if directory and not os.path.exists(directory):
            try:
                os.makedirs(directory, exist_ok=True)
            except OSError as e:
                raise FileOperationError(f"Cannot create directory {directory}: {e}")
        
        if os.path.exists(file_path) and not os.access(file_path, os.W_OK):
            raise FileOperationError(f"File not writable: {file_path}")
    
    @staticmethod
    def validate_file_size(file_path: str, max_size_mb: int = 100) -> None:
        """Validate that a file is not too large."""
        if os.path.exists(file_path):
            size_mb = os.path.getsize(file_path) / (1024 * 1024)
            if size_mb > max_size_mb:
                raise ValidationError(f"File too large: {size_mb:.1f}MB > {max_size_mb}MB")
    
    @staticmethod
    def validate_column_index(column_index: int, max_columns: int) -> None:
        """Validate column index."""
        if column_index < 0:
            raise ValidationError(f"Column index must be non-negative, got {column_index}")
        if column_index >= max_columns:
            raise ValidationError(f"Column index {column_index} out of range (max: {max_columns-1})")
    
    @staticmethod
    def validate_id_list(ids: List[str], id_type: str = "ID") -> None:
        """Validate a list of IDs."""
        if not ids:
            raise ValidationError(f"No {id_type}s provided")
        
        for i, id_val in enumerate(ids):
            if not id_val or not id_val.strip():
                raise ValidationError(f"Empty {id_type} at position {i}")
            if len(id_val.strip()) < 2:
                raise ValidationError(f"{id_type} too short: {id_val}")
    
    @staticmethod
    def validate_delimiter(delimiter: str) -> None:
        """Validate delimiter character."""
        if not delimiter:
            raise ValidationError("Delimiter cannot be empty")
        if len(delimiter) > 1:
            raise ValidationError(f"Delimiter must be single character, got: {delimiter}")
    
    @staticmethod
    def validate_encoding(encoding: str) -> None:
        """Validate encoding string."""
        try:
            "test".encode(encoding)
        except LookupError:
            raise ValidationError(f"Invalid encoding: {encoding}")
    
    @staticmethod
    def validate_plot_type(plot_type: str) -> None:
        """Validate plot type."""
        valid_types = ["bar", "hist", "scatter", "line"]
        if plot_type not in valid_types:
            raise ValidationError(f"Invalid plot type: {plot_type}. Valid types: {valid_types}")
    
    @staticmethod
    def validate_output_path(output_path: str, create_dir: bool = True) -> None:
        """Validate output path."""
        if not output_path:
            raise ValidationError("Output path cannot be empty")
        
        if create_dir:
            directory = os.path.dirname(output_path)
            if directory and not os.path.exists(directory):
                try:
                    os.makedirs(directory, exist_ok=True)
                except OSError as e:
                    raise FileOperationError(f"Cannot create output directory {directory}: {e}")
        
        DataValidator.validate_file_writable(output_path)
