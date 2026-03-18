
#!/usr/bin/env python3

import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_charge(filename):

    with h5py.File(filename, "r") as f:
        hits = f["Hits"]
        charge = hits["charge"]

        charge = charge[charge > 0]

        print("Entries:", len(charge))
        print("Mean:", np.mean(charge))
        print("Median:", np.median(charge))
        print("Std :", np.std(charge))

        # ヒスト
        plt.figure()
        plt.hist(charge, bins=200, range=(0, np.percentile(charge, 99.5)))
        plt.xlabel("Charge")
        plt.ylabel("Entries")
        plt.title("Charge distribution")
        plt.yscale("log")
        plt.grid()

        plt.savefig("charge_hist.png", dpi=150)
        print("Saved: charge_hist.png")

        # zoom版
        plt.figure()
        plt.hist(charge, bins=200, range=(0, 200))
        plt.xlabel("Charge")
        plt.ylabel("Entries")
        plt.title("Charge (zoom)")
        plt.yscale("log")
        plt.grid()

        plt.savefig("charge_hist_zoom.png", dpi=150)
        print("Saved: charge_hist_zoom.png")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python plot_charge.py file.h5")
        sys.exit(1)

    plot_charge(sys.argv[1])
