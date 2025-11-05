"""
Pytest configuration and fixtures for BETA tests
"""

import os
import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    """Return path to test data directory"""
    return Path(__file__).parent.parent / "BETA_test_data"


@pytest.fixture
def sample_peaks_file(test_data_dir):
    """Return path to sample peaks BED file"""
    return test_data_dir / "3656_peaks.bed"


@pytest.fixture
def sample_expr_file(test_data_dir):
    """Return path to sample expression file"""
    return test_data_dir / "ESR1_diff_expr.xls"


@pytest.fixture
def temp_output_dir(tmp_path):
    """Create temporary output directory"""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return output_dir
