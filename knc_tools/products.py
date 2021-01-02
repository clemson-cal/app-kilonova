import numpy as np
import pickle


primitive_dtype = np.dtype([('rho', float), ('ur', float), ('uq', float), ('pre', float)])


def fields():
    return list(primitive_dtype.fields) + ['entropy', 'scalar']


def array(d, dtype=float):
    return np.array(d['data'], dtype=dtype).reshape(d['dim'])


class Block:
    def __init__(self, radial_vertices=None, polar_vertices=None, scalar=None, primitive=None):
        self.radial_vertices = array(radial_vertices)
        self.polar_vertices  = array(polar_vertices)
        self.scalar          = array(scalar)
        self.primitive       = array(primitive, dtype=primitive_dtype)

    def field(self, key):
        if key in primitive_dtype.fields:
            return self.primitive[key]
        elif key == 'entropy':
            return self.field('pre') / self.field('rho')**(4/3)
        elif key == 'scalar':
            return self.scalar

    def pcolormesh_data(self, field, log=False):
        R, Q = [x.T for x in np.meshgrid(self.radial_vertices, self.polar_vertices)]
        x = R * np.sin(Q)
        z = R * np.cos(Q)
        c = self.field(field)
        return x, z, np.log10(c) if log else c


class Products:
    def __init__(self, filename):
        self._products = pickle.load(open(filename, 'rb'))

    @property
    def blocks(self):
        return self._products['blocks']

    @property
    def config(self):
        return self._products['config']

    @property
    def time(self):
        return self._products['time']
    
    def pcolormesh_data(self, field, log=False):
        return [Block(**b).pcolormesh_data(field, log=log) for b in self.blocks.values()]
