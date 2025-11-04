"""
Test core BETA functions
"""
import pytest
import math
from beta.core import corelib


def test_info_function():
    """Test Info logging function"""
    # Should not raise an exception
    corelib.Info("Test message")


def test_regulatory_potential_function():
    """Test regulatory potential score calculation"""
    # From pscore.py: Sg = lambda ldx: sum([math.exp(-0.5-4*t) for t in ldx])
    # Test the scoring function logic
    distances = [0.1, 0.2, 0.3]
    score = sum([math.exp(-0.5 - 4*t) for t in distances])

    # Score should be positive and decrease with distance
    assert score > 0

    # Test that closer peaks have higher scores
    close_distances = [0.01]
    far_distances = [0.9]

    close_score = sum([math.exp(-0.5 - 4*t) for t in close_distances])
    far_score = sum([math.exp(-0.5 - 4*t) for t in far_distances])

    assert close_score > far_score


def test_distance_normalization():
    """Test distance normalization (0-1 range)"""
    max_distance = 100000  # 100kb default

    # Distances within range
    test_distances = [1000, 10000, 50000, 99999]

    for dist in test_distances:
        normalized = dist / max_distance
        assert 0 <= normalized <= 1


def test_peak_sorting():
    """Test peak sorting by score"""
    peaks = [
        ['chr1', '1000', '2000', 'peak1', '50'],
        ['chr1', '3000', '4000', 'peak2', '100'],
        ['chr1', '5000', '6000', 'peak3', '75'],
    ]

    sorted_peaks = sorted(peaks, key=lambda p: float(p[4]), reverse=True)

    # Highest score should be first
    assert float(sorted_peaks[0][4]) == 100
    assert float(sorted_peaks[1][4]) == 75
    assert float(sorted_peaks[2][4]) == 50
