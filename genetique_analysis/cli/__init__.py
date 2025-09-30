"""
Command-line interface for genetique analysis.

This module provides the command-line interface for running genetics analysis
from the terminal.
"""

from .error_tests import cli_error_tests, main_error_recaptures
from .main import load_configuration, main

__all__ = ["main", "load_configuration", "main_error_recaptures", "cli_error_tests"]
