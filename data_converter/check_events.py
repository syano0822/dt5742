#!/usr/bin/env python3

import h5py
import numpy as np
import sys


def check_events(filename):
    with h5py.File(filename, "r") as f:
        hits = f["Hits"]

        if "event" not in hits.dtype.names:
            print("ERROR: 'event' field not found")
            return

        events = hits["event"]

        unique, counts = np.unique(events, return_counts=True)

        print("=== Event statistics ===")
        print("Number of events:", len(unique))
        print("Mean hits/event:", np.mean(counts))
        print("Max hits/event:", np.max(counts))
        print("Min hits/event:", np.min(counts))

        print("\n=== First 10 events ===")
        for i in range(min(10, len(unique))):
            print(f"event {unique[i]} : {counts[i]} hits")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python check_events.py file.h5")
        sys.exit(1)

    check_events(sys.argv[1])
