"""
Utility functions for BETA
Modernized version with proper logging support
"""

import sys
import os
import time
import logging
import subprocess
from subprocess import call as subpcall
from typing import Optional


# Setup logging
def setup_logger(name: str = "beta", level: int = logging.INFO) -> logging.Logger:
    """
    Setup a logger with consistent formatting

    Args:
        name: Logger name
        level: Logging level

    Returns:
        Configured logger instance
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Remove existing handlers
    logger.handlers.clear()

    # Create console handler with formatting
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(level)

    # Create formatter
    formatter = logging.Formatter(
        '[%(asctime)s] %(levelname)s: %(message)s',
        datefmt='%H:%M:%S'
    )
    handler.setFormatter(formatter)

    logger.addHandler(handler)
    return logger


# Create default logger
_logger = setup_logger()

# Function aliases for backward compatibility
error = _logger.critical
warn = _logger.warning
info = _logger.info
debug = _logger.debug


def Info(info_str: str) -> None:
    """
    Print information with timestamp (backward compatible function)

    Args:
        info_str: Information string to print
    """
    print(f"[{time.strftime('%H:%M:%S')}] {info_str}")


def run_cmd(command: str) -> int:
    """
    Execute a shell command

    Args:
        command: Shell command to execute

    Returns:
        Return code from command execution
    """
    return subpcall(command, shell=True)


def run_pip(command: str) -> subprocess.Popen:
    """
    Execute a command and return Popen object for capturing output

    Args:
        command: Shell command to execute

    Returns:
        Popen object for the command
    """
    return subprocess.Popen(
        command,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def format_time(seconds: float) -> str:
    """
    Format seconds into HH:MM:SS format

    Args:
        seconds: Time in seconds

    Returns:
        Formatted time string
    """
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)
    return f"{hours}:{minutes:02d}:{secs:02d}"


def print_elapsed_time(start_time: float) -> None:
    """
    Print elapsed time since start_time

    Args:
        start_time: Start time from time.time()
    """
    elapsed = time.time() - start_time
    formatted_time = format_time(elapsed)
    Info(f"Total time: {formatted_time}")
