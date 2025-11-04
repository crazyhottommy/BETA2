# PyPI Release Guide for BETA

This guide walks you through publishing BETA to PyPI (Python Package Index).

## Prerequisites

### 1. Create PyPI Accounts
You'll need accounts on both Test PyPI and Production PyPI:

- **Test PyPI** (for testing): https://test.pypi.org/account/register/
- **Production PyPI** (for release): https://pypi.org/account/register/

### 2. Install Required Tools
```bash
pip install build twine
```

### 3. Setup API Tokens (Recommended)

#### Test PyPI Token
1. Go to https://test.pypi.org/manage/account/token/
2. Create a new API token with scope "Entire account"
3. Copy the token (starts with `pypi-`)
4. Save it securely

#### Production PyPI Token
1. Go to https://pypi.org/manage/account/token/
2. Create a new API token with scope "Entire account"
3. Copy the token
4. Save it securely

### 4. Configure PyPI Credentials

Create `~/.pypirc` file:
```ini
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
username = __token__
password = pypi-YOUR_PRODUCTION_TOKEN_HERE

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-YOUR_TEST_TOKEN_HERE
```

**Important**: Keep this file secure! Consider using `chmod 600 ~/.pypirc`

## Pre-Release Checklist

Before publishing, verify:

- [ ] All tests pass: `pytest tests/ -v`
- [ ] Package builds: `python -m build`
- [ ] Package checks pass: `twine check dist/*`
- [ ] Version number is correct in `pyproject.toml`
- [ ] README.md is up to date
- [ ] CHANGELOG.md exists and is current
- [ ] LICENSE file exists
- [ ] GitHub repository URL is correct in `pyproject.toml`
- [ ] All changes are committed to git
- [ ] Git tag created: `git tag v2.0.0`

## Release Process

### Step 1: Clean Build Environment
```bash
# Remove old builds
rm -rf dist/ build/ *.egg-info

# Ensure clean state
git status
```

### Step 2: Build the Package
```bash
python -m build
```

This creates:
- `dist/beta_binding_analysis-2.0.0-py3-none-any.whl` (wheel)
- `dist/beta_binding_analysis-2.0.0.tar.gz` (source distribution)

### Step 3: Verify the Build
```bash
# Check package integrity
twine check dist/*

# Inspect contents
tar -tzf dist/beta_binding_analysis-2.0.0.tar.gz | head -20
```

### Step 4: Test on Test PyPI (RECOMMENDED)

Upload to Test PyPI first to catch any issues:
```bash
twine upload --repository testpypi dist/*
```

After upload, test the installation:
```bash
# Create a fresh virtual environment
python -m venv test_env
source test_env/bin/activate  # On Windows: test_env\Scripts\activate

# Install from Test PyPI
pip install --index-url https://test.pypi.org/simple/ \
    --extra-index-url https://pypi.org/simple/ \
    beta-binding-analysis

# Test it works
beta --version
beta --help

# Cleanup
deactivate
rm -rf test_env
```

### Step 5: Upload to Production PyPI

If Test PyPI installation works perfectly:
```bash
twine upload dist/*
```

You'll be prompted for your username and password (use `__token__` and your API token).

### Step 6: Verify Production Installation

Test the production installation:
```bash
# Create a fresh virtual environment
python -m venv prod_test
source prod_test/bin/activate

# Install from PyPI
pip install beta-binding-analysis

# Test it
beta --version
beta --help

# Cleanup
deactivate
rm -rf prod_test
```

### Step 7: Create GitHub Release

1. Go to your GitHub repository
2. Click "Releases" → "Create a new release"
3. Tag: `v2.0.0`
4. Title: `BETA v2.0.0 - Python 3 Release`
5. Description: Copy from CHANGELOG.md
6. Attach the wheel and tar.gz files from `dist/`
7. Publish release

## Post-Release

### Update Version for Development
In `pyproject.toml`, bump version to next development version:
```toml
version = "2.0.1.dev0"
```

### Announce the Release
- Update README badges
- Post on relevant forums/mailing lists
- Update documentation site
- Tweet/social media announcement

## Troubleshooting

### Issue: "File already exists"
If you get this error, you've already uploaded this version. You must:
1. Increment the version number in `pyproject.toml`
2. Rebuild: `rm -rf dist/ && python -m build`
3. Upload again

### Issue: Package not found when installing
Wait a few minutes. PyPI can take time to index new packages.

### Issue: Import errors after installation
Check that:
- All dependencies are in `pyproject.toml`
- Package data files are included
- C extensions compiled correctly

### Issue: README not rendering on PyPI
PyPI only renders certain Markdown features. Check:
- Use standard CommonMark Markdown
- Validate with: `twine check dist/*`
- Preview at: https://github.com/pypa/readme_renderer

## Package Maintenance

### Releasing Updates

For patch releases (2.0.1, 2.0.2):
1. Update version in `pyproject.toml`
2. Update CHANGELOG.md
3. Follow steps 1-7 above

For minor releases (2.1.0, 2.2.0):
1. Update version
2. Update CHANGELOG.md
3. May add new features
4. Follow steps 1-7

For major releases (3.0.0):
1. Update version
2. Update CHANGELOG.md
3. May include breaking changes
4. Follow steps 1-7
5. Announce breaking changes prominently

### Yanking a Release

If you need to remove a broken release from PyPI:
```bash
# This doesn't delete, but marks it as unavailable
twine upload --skip-existing dist/*
# Then mark the bad version as "yanked" in PyPI web interface
```

## Security Considerations

1. **Never commit PyPI tokens to git**
2. **Use API tokens instead of passwords**
3. **Enable 2FA on PyPI account**
4. **Regularly rotate API tokens**
5. **Use separate tokens for Test and Production PyPI**
6. **Set token scope to specific projects when possible**

## Useful Commands

```bash
# Check if package name is available
pip install beta-binding-analysis  # If it fails, name is available

# Check package on PyPI
pip index versions beta-binding-analysis

# Download package without installing
pip download beta-binding-analysis --no-deps

# View package metadata
pip show beta-binding-analysis

# List all versions
pip index versions beta-binding-analysis
```

## References

- PyPI: https://pypi.org/
- Test PyPI: https://test.pypi.org/
- Packaging Guide: https://packaging.python.org/
- Twine: https://twine.readthedocs.io/
- Build: https://build.pypa.io/

## Quick Reference Card

```bash
# Complete release workflow
rm -rf dist/ build/ *.egg-info     # Clean
python -m build                      # Build
twine check dist/*                   # Verify
twine upload --repository testpypi dist/*  # Test
twine upload dist/*                  # Release
git tag v2.0.0 && git push --tags   # Tag
```

---

**Package Name**: `beta-binding-analysis`
**PyPI URL**: https://pypi.org/project/beta-binding-analysis/
**Documentation**: http://cistrome.org/BETA/tutorial.html
