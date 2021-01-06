import numpy as np
import pickle


primitive_dtype = np.dtype([('rho', float), ('ur', float), ('uq', float), ('pre', float)])


def fields():
    return list(primitive_dtype.fields) + ['entropy', 'scalar']


def array(d, dtype=float):
    if type(d['data'][0]) is float:
        data = d['data']
    elif type(d['data'][0]) is tuple:
        data = d['data']
    elif type(d['data'][0] is list):
        data = [tuple(t) for t in d['data']]
    else:
        raise IOError('unknown array element format')
    return np.array(data, dtype=dtype).reshape(d['dim'])


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
        file = open(filename, 'rb')

        if filename.endswith('.pk'):
            self._products = pickle.load(file)
        elif filename.endswith('.cbor'):
            import cbor2            
            self._products = cbor2.load(file)
        elif filename.endswith('.cboz'):
            import snappy, cbor2
            stream = snappy.StreamDecompressor()
            data = stream.decompress(file.read())
            self._products = cbor2.loads(data)
        else:
            raise IOError(f'unknown file format {filename}')

    @property
    def blocks(self):
        return [Block(**self._products['blocks'][k]) for k in sorted(self._products['blocks'].keys())]

    @property
    def config(self):
        return self._products['config']

    @property
    def time(self):
        return self._products['time']
    
    def pcolormesh_data(self, field, log=False):
        return [b.pcolormesh_data(field, log=log) for b in self.blocks]

    def radial_profile(self, field, polar_index=0):
        return np.concatenate([b.field(field)[:, polar_index] for b in self.blocks])

    def radial_vertices(self):
        return np.concatenate([b.radial_vertices[:-1] for b in self.blocks])
