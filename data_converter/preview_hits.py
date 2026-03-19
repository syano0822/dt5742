#!/usr/bin/env python3

import h5py
import numpy as np
import sys


def preview_hits(filename, max_events=5):
    with h5py.File(filename, "r") as f:

        if "Hits" not in f:
            print("ERROR: No 'Hits' dataset found")
            return

        hits = f["Hits"]

        print("=== Dataset info ===")
        print("shape:", hits.shape)
        print("dtype:", hits.dtype)

        print("\n=== First entries ===")

        for i in range(min(max_events, len(hits))):
            row = hits[i]
            print(f"{i}: {row}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python preview_hits.py file.h5")
        sys.exit(1)

    preview_hits(sys.argv[1])
