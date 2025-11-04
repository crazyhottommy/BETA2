# BETA Test Suite

This directory contains the test suite for BETA (Binding and Expression Target Analysis).

## Test Structure

- `test_cli.py` - Tests for command-line interface
- `test_imports.py` - Tests for module imports
- `test_core_functions.py` - Tests for core functionality (scoring, distance calculations)
- `test_fileformat_check.py` - Tests for file format validation
- `conftest.py` - Pytest fixtures and configuration

## Running Tests

### Run all tests
```bash
pytest tests/
```

### Run with coverage
```bash
pytest --cov=beta --cov-report=html tests/
```

### Run specific test file
```bash
pytest tests/test_cli.py
```

### Run specific test
```bash
pytest tests/test_cli.py::test_beta_version -v
```

## Requirements

The test suite requires:
- pytest >= 7.0
- pytest-cov >= 3.0

Install with:
```bash
pip install -e ".[dev]"
```

## Test Data

Tests use data from the `BETA_test_data/` directory, which includes:
- Sample ChIP-seq peak files
- Sample expression files
- Reference genome annotations
