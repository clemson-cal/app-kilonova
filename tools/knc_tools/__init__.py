from .knc_tools import app, products
from .timed import timed


def load_products(filename):
    if filename.startswith('prods'):
        with timed('load products'):
            p = products(filename)
    else:
        with timed('load app'):
            a = app(filename)
        with timed('make products'):
            p = a.make_products()
    return p
