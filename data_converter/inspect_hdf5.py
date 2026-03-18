#!/usr/bin/env python3

import h5py
import sys


def print_structure(name, obj):
    print(f"{name} ({type(obj).__name__})")
    if isinstance(obj, h5py.Dataset):
        print(f"  shape: {obj.shape}")
        print(f"  dtype: {obj.dtype}")
        if obj.attrs:
            print("  attrs:")
            for k, v in obj.attrs.items():
                print(f"    {k}: {v}")


def inspect_hdf5(filename):
    print("=" * 60)
    print(f"File: {filename}")
    print("=" * 60)

    with h5py.File(filename, "r") as f:
        print("\n--- Structure ---")
        f.visititems(print_structure)

        print("\n--- Top-level keys ---")
        for key in f.keys():
            print(f"  {key}")

        print("\n--- Attributes ---")
        for k, v in f.attrs.items():
            print(f"  {k}: {v}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python inspect_hdf5.py file.h5")
        sys.exit(1)

    inspect_hdf5(sys.argv[1])
