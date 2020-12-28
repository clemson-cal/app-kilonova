import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt


primitive_dtype = np.dtype([('rho', float), ('ur', float), ('uq', float), ('pre', float)])


class Block:
    def __init__(self, radial_vertices=None, polar_vertices=None, scalar=None, primitive=None):
        self.radial_vertices = array(radial_vertices)
        self.polar_vertices  = array(polar_vertices)
        self.scalar          = array(scalar)
        self.primitive       = array(primitive, dtype=primitive_dtype)

    def pcolormesh_data(self, field, log=False):
        R, Q = [x.T for x in np.meshgrid(self.radial_vertices, self.polar_vertices)]
        x = R * np.sin(Q)
        z = R * np.cos(Q)
        c = self.primitive[field]
        return x, z, np.log10(c) if log else c


def array(d, dtype=float):
    return np.array(d['data'], dtype=dtype).reshape(d['dim'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("filenames", nargs='+')
    parser.add_argument('-r', '--range', default='None,None', help='vmin and vmax parameters for the relief plot')
    parser.add_argument('-f', '--field', default='rho', choices=primitive_dtype.fields)
    parser.add_argument('-l', '--log', action='store_true')
    args = parser.parse_args()

    vmin, vmax = eval(args.range)

    for filename in args.filenames:
        fig = plt.figure(figsize=[8, 12])
        ax1 = fig.add_subplot(1, 1, 1)
        prods = pickle.load(open(filename, 'rb'))

        pcolormesh_data = [Block(**block).pcolormesh_data(args.field, log=args.log) for block in prods['blocks'].values()]
        vmin = min([c.min() for _, _, c in pcolormesh_data]) if vmin is None else vmin
        vmax = max([c.max() for _, _, c in pcolormesh_data]) if vmax is None else vmax

        for x, z, c in pcolormesh_data:
            cm = ax1.pcolormesh(x, z, c, vmin=vmin, vmax=vmax, edgecolors='none')

        ax1.set_aspect('equal')
        fig.colorbar(cm)
    plt.show()


if __name__ == "__main__":
    main()
