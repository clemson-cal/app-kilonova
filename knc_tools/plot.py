import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import knc_loader
from timed import timed


def load_products(filename):
    if 'prods' in filename:
        with timed('load products'):
            p = knc_loader.products(filename)
    else:
        with timed('load app'):
            a = knc_loader.app(filename)
        with timed('make products'):
            p = a.make_products()
        with timed('cache products'):
            p.save(filename.replace('chkpt', 'prods'))
    return p


def block_vertices(block):
    R, Q = [x.T for x in np.meshgrid(block.radial_vertices, block.polar_vertices)]
    return R * np.sin(Q), R * np.cos(Q)


def mesh_vertices(products):
    return [block_vertices(products[k]) for k in products]


def known_fields():
    ['rho', 'pre', 'ur', 'uq', 'scalar']


def block_field(block, field, transform=lambda x: x):
    if field == 'rho':
        d = block.comoving_mass_density
    elif field == 'pre':
        d = block.gas_pressure
    elif field == 'ur':
        d = block.radial_four_velocity
    elif field == 'uq':
        d = block.polar_four_velocity
    elif field == 'scalar':
        d = block.scalar
    else:
        raise ValueError(f'unknown field {field}')
    return transform(d)


def mesh_field(products, field, transform=lambda x: x):
    return [block_field(products[k], field, transform) for k in products]


def variable(args):
    if args.field == 'rho' and args.log:
        return r'$\log_{{10}}(\rho) \ [\rm{{g/cm^3}}]$'
    else:
        return args.field


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    parser.add_argument('-r', '--range', default='None,None', help='vmin and vmax parameters for the relief plot')
    parser.add_argument('-f', '--field', default='rho', choices=known_fields())
    parser.add_argument('-l', '--log', action='store_true')
    parser.add_argument('-c', '--cmap', default='viridis')
    parser.add_argument('--hardcopy', action='store_true')
    parser.add_argument('--radius', default=None, type=float)
    args = parser.parse_args()

    vmin, vmax = eval(args.range)
    fig = plt.figure(figsize=[7, 11])
    ax1 = fig.add_subplot(1, 1, 1)

    filename = args.filenames[0]
    products = load_products(filename)

    with timed('make vertices and fields'):
        vertices = mesh_vertices(products)
        field = mesh_field(products, args.field, np.log10 if args.log else lambda x: x)
        setup = next(iter(products.config['model']))
        vmin = min([c.min() for c in field])
        vmax = max([c.max() for c in field])

    with timed('load plots'):
        for (x, z), c in zip(vertices, field):
            cm = ax1.pcolormesh(x, z, c, vmin=vmin, vmax=vmax, edgecolors='none', cmap=args.cmap)

    if args.radius is not None:
        ax1.set_xlim(0, args.radius)
        ax1.set_ylim(-args.radius, args.radius)
    ax1.set_xlabel(r'$x \ [\rm{cm}]$')
    ax1.set_ylabel(r'$z \ [\rm{cm}]$')

    ax1.set_aspect('equal')
    fig.colorbar(cm)
    fig.subplots_adjust(left=0, right=1, top=0.9, bottom=0.05)
    fig.suptitle(r'Setup: $\mathtt{{{}}}$   {}   $t = {:.4}s$'.format(setup.replace('_', '-'), variable(args), products.time))

    with timed('show'):
        plt.show()


if __name__ == "__main__":
    main()
