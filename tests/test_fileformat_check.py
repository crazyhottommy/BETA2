"""
Test file format checking functionality
"""
import pytest
from pathlib import Path


def test_peak_file_exists(sample_peaks_file):
    """Test that sample peak file exists"""
    assert sample_peaks_file.exists()


def test_peak_file_format(sample_peaks_file):
    """Test that peak file has correct BED format"""
    with open(sample_peaks_file, 'r') as f:
        first_line = f.readline().strip()
        fields = first_line.split('\t')

        # Should have at least 3 columns (chr, start, end)
        assert len(fields) >= 3

        # First column should be chromosome
        assert fields[0].startswith('chr')

        # Second and third columns should be numeric (positions)
        assert fields[1].isdigit()
        assert fields[2].isdigit()

        # Start should be less than end
        assert int(fields[1]) < int(fields[2])


def test_expr_file_exists(sample_expr_file):
    """Test that sample expression file exists"""
    assert sample_expr_file.exists()


def test_expr_file_has_header(sample_expr_file):
    """Test that expression file has proper header"""
    with open(sample_expr_file, 'r') as f:
        header = f.readline().strip()
        # Should have tab-separated columns
        assert '\t' in header
        fields = header.split('\t')
        assert len(fields) >= 3
