import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
import products


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    parser.add_argument('-r', '--range', default='None,None', help='vmin and vmax parameters for the relief plot')
    parser.add_argument('-f', '--field', default='rho', choices=products.fields())
    parser.add_argument('-l', '--log', action='store_true')
    parser.add_argument('--radius', default=None, type=float)
    args = parser.parse_args()

    vmin, vmax = eval(args.range)

    for filename in args.filenames:
        fig = plt.figure(figsize=[7, 10])
        ax1 = fig.add_subplot(1, 1, 1)

        prods = products.Products(filename)
        pcolormesh_data = prods.pcolormesh_data(args.field, log=args.log)
        vmin = min([c.min() for _, _, c in pcolormesh_data]) if vmin is None else vmin
        vmax = max([c.max() for _, _, c in pcolormesh_data]) if vmax is None else vmax

        for x, z, c in pcolormesh_data:
            cm = ax1.pcolormesh(x, z, c, vmin=vmin, vmax=vmax, edgecolors='none')

        if args.radius is not None:
            ax1.set_xlim(0, args.radius)
            ax1.set_ylim(-args.radius, args.radius)
        ax1.set_xlabel(r'$x \ [\rm{cm}]$')
        ax1.set_ylabel(r'$z \ [\rm{cm}]$')
        ax1.set_aspect('equal')
        fig.colorbar(cm)
        fig.suptitle(r'$t = {:.4}s$'.format(prods.time))
    plt.show()


if __name__ == "__main__":
    main()
