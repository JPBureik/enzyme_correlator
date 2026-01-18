#!/usr/bin/env python3
"""Backward compatibility wrapper for enzyme_correlator.

This file is kept for backward compatibility with the original invocation method.
For new usage, install the package and run: enzyme-correlator
"""

from enzyme_correlator import main

if __name__ == "__main__":
    main()
