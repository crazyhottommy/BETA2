# GitHub Actions Workflows

This directory contains CI/CD workflows for BETA.

## Workflows

### Tests (`tests.yml`)
- **Triggers**: Push and PR to `main` or `develop`
- **Matrix**: Python 3.8-3.12 on Ubuntu and macOS
- **Actions**:
  - Install dependencies
  - Compile MISP C code
  - Run pytest with coverage
  - Upload coverage to Codecov

### Build (`build.yml`)
- **Triggers**: Push to `main`, PRs, releases
- **Actions**:
  - Build Python package (wheel and sdist)
  - Check package with twine
  - Test installation on multiple platforms

### Lint (`lint.yml`)
- **Triggers**: Push and PR to `main` or `develop`
- **Actions**:
  - Check code formatting with black
  - Lint with flake8
  - Type checking with mypy (non-blocking)

## Badges

Add these to your README.md:

```markdown
[![Tests](https://github.com/yourusername/BETA2/actions/workflows/tests.yml/badge.svg)](https://github.com/yourusername/BETA2/actions/workflows/tests.yml)
[![Build](https://github.com/yourusername/BETA2/actions/workflows/build.yml/badge.svg)](https://github.com/yourusername/BETA2/actions/workflows/build.yml)
[![Lint](https://github.com/yourusername/BETA2/actions/workflows/lint.yml/badge.svg)](https://github.com/yourusername/BETA2/actions/workflows/lint.yml)
```

## Local Testing

Test the workflows locally before pushing:

```bash
# Run tests
pytest tests/ -v --cov=beta

# Check formatting
black --check src/beta tests/

# Lint code
flake8 src/beta --max-line-length=100

# Build package
python -m build
twine check dist/*
```
