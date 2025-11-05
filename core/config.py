import os
import yaml
from dataclasses import dataclass, field
from typing import Dict, Any, Optional
from pathlib import Path


@dataclass(frozen=True)
class AppConfig:
    app_name: str = "LeishDomains"
    version: str = "2.0.0"
    default_session_prefix: str = "SESSION"
    config_file: Optional[str] = None
    _config_data: Dict[str, Any] = field(default_factory=dict, init=False)

    def __post_init__(self):
        """Load configuration from file if specified."""
        if self.config_file and os.path.exists(self.config_file):
            try:
                with open(self.config_file, 'r', encoding='utf-8') as f:
                    config_data = yaml.safe_load(f)
                    object.__setattr__(self, '_config_data', config_data or {})
            except Exception:
                # Fall back to defaults if config file can't be loaded
                object.__setattr__(self, '_config_data', {})

    @property
    def cwd(self) -> str:
        return os.getcwd()

    @property
    def log_level(self) -> str:
        return self._config_data.get('logging', {}).get('level', 'INFO')

    @property
    def log_format(self) -> str:
        return self._config_data.get('logging', {}).get('format', 
            '%(asctime)s | %(levelname)s | %(name)s | %(message)s')

    @property
    def log_file(self) -> str:
        return self._config_data.get('logging', {}).get('file', 'leishdomains.log')

    @property
    def max_file_size_mb(self) -> int:
        return self._config_data.get('data', {}).get('max_file_size_mb', 100)

    @property
    def default_delimiter(self) -> str:
        return self._config_data.get('data', {}).get('default_delimiter', '\t')

    @property
    def encoding(self) -> str:
        return self._config_data.get('data', {}).get('encoding', 'utf-8')

    @classmethod
    def from_file(cls, config_path: str) -> 'AppConfig':
        """Create AppConfig from configuration file."""
        return cls(config_file=config_path)

    @classmethod
    def default(cls) -> 'AppConfig':
        """Create AppConfig with default values."""
        default_config = Path(__file__).parent.parent / "config" / "default.yaml"
        if default_config.exists():
            return cls.from_file(str(default_config))
        return cls()


