import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
import products


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    args = parser.parse_args()

    p = products.Products(args.filename)
    r = p.radial_vertices()
    s = p.radial_profile('scalar', polar_index=0)
    u = p.radial_profile('ur', polar_index=0)
    i = np.where(u > 1)[0][0]

    fig = plt.figure(figsize=[8, 10])
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    ax1.plot(r, s)
    ax2.plot(r, u)
    ax1.axhline(s[i], lw=0.5, c='k')
    ax2.axhline(u[i], lw=0.5, c='k')
    ax1.axvline(r[i], lw=0.5, c='k')
    ax2.axvline(r[i], lw=0.5, c='k')
    ax1.set_ylabel(r'Mass coordinate [$\rm{g}$]')
    ax2.set_ylabel(r'Radial $\Gamma \beta$')
    ax2.set_xlabel(r'$r \ [\rm{cm}]$')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    plt.show()


if __name__ == "__main__":
    main()
