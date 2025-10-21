"""Statistical analysis and data processing."""

import numpy as np
from scipy import stats
from collections import Counter
from typing import List, Dict, Tuple, Optional, Any
import decimal


class DataAnalyzer:
    """Handles statistical analysis of data."""
    
    def __init__(self):
        self.logger = None  # Will be set by caller if needed
    
    def _log(self, message: str):
        """Log message if logger is available."""
        if self.logger:
            self.logger.info(message)
        else:
            print(message)
    
    def count_frequencies(self, data: List[Any]) -> Dict[Any, int]:
        """Count frequency of items in data."""
        return dict(Counter(data))
    
    def find_most_common(self, data: List[Any], n: int = 10) -> List[Tuple[Any, int]]:
        """Find the n most common items in data."""
        counter = Counter(data)
        return counter.most_common(n)
    
    def calculate_statistics(self, values: List[float]) -> Dict[str, float]:
        """Calculate basic statistics for a list of values."""
        if not values:
            return {}
        
        values_array = np.array(values)
        
        stats_dict = {
            'count': len(values),
            'mean': np.mean(values_array),
            'median': np.median(values_array),
            'std': np.std(values_array),
            'min': np.min(values_array),
            'max': np.max(values_array),
            'sum': np.sum(values_array)
        }
        
        return stats_dict
    
    def remove_exponent(self, value: Any) -> str:
        """Remove scientific notation from decimal values."""
        decimal_places = 8
        max_digits = 16
        
        if isinstance(value, decimal.Decimal):
            context = decimal.getcontext().copy()
            context.prec = max_digits
            return "{0:f}".format(value.quantize(decimal.Decimal(".1") ** decimal_places, context=context))
        else:
            return "%.*f" % (decimal_places, value)
    
    def is_number(self, s: str) -> bool:
        """Check if a string represents a number."""
        try:
            float(s)
            return True
        except ValueError:
            return False
    
    def analyze_column(self, filename: str, column_index: int, 
                      delimiter: str = '\t') -> Dict[str, Any]:
        """Analyze a specific column in a file."""
        values = []
        non_numeric_count = 0
        
        with open(filename, 'r', encoding='utf-8') as f:
            # Skip header
            f.readline()
            
            for line in f:
                line = line.strip()
                if line:
                    items = line.split(delimiter)
                    if len(items) > column_index:
                        value = items[column_index].strip()
                        if self.is_number(value):
                            values.append(float(value))
                        else:
                            non_numeric_count += 1
        
        if not values:
            return {'error': 'No numeric values found'}
        
        stats = self.calculate_statistics(values)
        stats['non_numeric_count'] = non_numeric_count
        
        return stats
