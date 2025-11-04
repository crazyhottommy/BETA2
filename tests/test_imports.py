"""
Test that all BETA modules can be imported without errors
"""
import pytest


def test_import_beta():
    """Test main beta package imports"""
    import beta
    assert hasattr(beta, '__version__')


def test_import_cli():
    """Test CLI module imports"""
    from beta import cli
    assert hasattr(cli, 'main')
    assert hasattr(cli, 'prepare_argparser')


def test_import_core_modules():
    """Test core module imports"""
    from beta.core import corelib
    from beta.core import pscore
    from beta.core import up_down_distance
    from beta.core import up_down_score
    from beta.core import expr_combine
    from beta.core import fileformat_check
    from beta.core import permp
    from beta.core import runbeta
    from beta.core import opt_validator

    assert hasattr(corelib, 'Info')
    assert hasattr(corelib, 'run_cmd')


def test_import_motif_modules():
    """Test motif module imports"""
    from beta.motif import motif_parser
    from beta.motif import motif_scan
    from beta.motif import motif_clustering
    from beta.motif import fastafrombed
    from beta.motif import bayesian_motif_comp


def test_import_utils():
    """Test utils module imports"""
    from beta import utils
