#!/usr/bin/env python3
"""
Setup script for BETA (Binding and Expression Target Analysis)
Maintained for backward compatibility. Modern installations should use pyproject.toml.
"""

import os
import sys
from pathlib import Path
from setuptools import setup, find_packages
from subprocess import call as subpcall

# Ensure Python 3.8+
if sys.version_info < (3, 8):
    sys.stderr.write("ERROR: Python version must be 3.8 or higher!\n")
    sys.stderr.write(f"Current version: {sys.version}\n")
    sys.exit(1)


def compile_misp():
    """Compile the MISP C program for motif scanning"""
    curdir = os.getcwd()
    misp_dir = Path('src/beta/misp')

    if misp_dir.exists():
        os.chdir(misp_dir)
        print("Compiling MISP motif scanner...")
        try:
            subpcall('make', shell=True)
            subpcall('chmod 755 misp', shell=True)
            print("MISP compiled successfully")
        except Exception as e:
            print(f"Warning: Could not compile MISP: {e}")
            print("You may need to compile it manually later")
        finally:
            os.chdir(curdir)


if __name__ == '__main__':
    # Compile MISP if needed
    if '--skip-misp' not in sys.argv:
        compile_misp()
    else:
        sys.argv.remove('--skip-misp')

    setup()
