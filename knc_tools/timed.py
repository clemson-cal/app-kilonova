import time


class timed(object): 
    def __init__(self, message):
        self.message = message
        self.start = time.perf_counter()

    def __enter__(self): 
        print('{:.<32} '.format(self.message), end='', flush=True)

    def __exit__(self, *args): 
        print('{:.3}s'.format(time.perf_counter() - self.start))
