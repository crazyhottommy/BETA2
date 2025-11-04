"""
BETA: Binding and Expression Target Analysis

A tool for integrative analysis of ChIP-seq and RNA-seq/microarray data
to predict transcription factor direct target genes.
"""

__version__ = "2.0.0"
__author__ = "Su Wang, Tommy Tang"
__email__ = "wangsu0623@gmail.com"

from beta.utils import setup_logger, Info, error, warn

__all__ = ['__version__', 'setup_logger', 'Info', 'error', 'warn']
