#!/usr/bin/env python3

import h5py
import sys


REQUIRED_FIELDS = ["event", "column", "row"]


def check_corry_format(filename):
    with h5py.File(filename, "r") as f:

        if "Hits" not in f:
            print("❌ Missing dataset: Hits")
            return

        hits = f["Hits"]

        print("Fields:", hits.dtype.names)

        for field in REQUIRED_FIELDS:
            if field not in hits.dtype.names:
                print(f"❌ Missing field: {field}")
            else:
                print(f"✔ Found: {field}")

        if "charge" in hits.dtype.names:
            print("✔ charge field present")
        else:
            print("⚠ charge not found (OK for binary detectors)")

        if "timestamp" in hits.dtype.names:
            print("✔ timestamp field present")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_corry_format.py file.h5")
        sys.exit(1)

    check_corry_format(sys.argv[1])
