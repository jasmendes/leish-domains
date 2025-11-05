"""Custom exceptions for LeishDomains."""


class LeishDomainsError(Exception):
    """Base exception for LeishDomains."""
    pass


class ConfigurationError(LeishDomainsError):
    """Raised when there's a configuration error."""
    pass


class SessionError(LeishDomainsError):
    """Raised when there's a session management error."""
    pass


class DataProcessingError(LeishDomainsError):
    """Raised when there's a data processing error."""
    pass


class FileOperationError(LeishDomainsError):
    """Raised when there's a file operation error."""
    pass


class AnalysisError(LeishDomainsError):
    """Raised when there's an analysis error."""
    pass


class ValidationError(LeishDomainsError):
    """Raised when data validation fails."""
    pass
