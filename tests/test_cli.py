"""
Test BETA command-line interface
"""
import pytest
import subprocess
import sys


def test_beta_version():
    """Test that beta --version works"""
    result = subprocess.run(
        ["beta", "--version"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "2.0.0" in result.stdout


def test_beta_help():
    """Test that beta --help works"""
    result = subprocess.run(
        ["beta", "--help"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "BETA" in result.stdout
    assert "basic" in result.stdout
    assert "plus" in result.stdout
    assert "minus" in result.stdout


def test_beta_basic_help():
    """Test that beta basic --help works"""
    result = subprocess.run(
        ["beta", "basic", "--help"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "peakfile" in result.stdout
    assert "diff_expr" in result.stdout


def test_beta_plus_help():
    """Test that beta plus --help works"""
    result = subprocess.run(
        ["beta", "plus", "--help"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "peakfile" in result.stdout
    assert "genome_sequence" in result.stdout or "gs" in result.stdout


def test_beta_minus_help():
    """Test that beta minus --help works"""
    result = subprocess.run(
        ["beta", "minus", "--help"],
        capture_output=True,
        text=True
    )
    assert result.returncode == 0
    assert "peakfile" in result.stdout
