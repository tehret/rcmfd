import os
import cffi
import numpy as np

root = os.path.dirname(os.path.realpath(__file__))
ffi = cffi.FFI()
with open(os.path.join(root, 'build/all.h')) as header:
    ffi.cdef(header.read())

rcmfd = ffi.dlopen(os.path.join(root, "build/libwrcmfd.so"))

def P(array):
    typestr = 'double*'
    if array.dtype == np.float32:
        typestr = 'float*'
    elif array.dtype == np.bool:
        typestr = 'bool*'
    elif array.dtype == np.int32:
        typestr = 'int*'
    # requires cffi 0.12
    return ffi.from_buffer(typestr, array, require_writable=True)

def STR(ptr):
    return ffi.string(ptr)

class a(object):
    pass

def wrap(a):
    if isinstance(a, np.ndarray):
        return P(a)
    elif isinstance(a, list):
        return list(map(wrap, a))
    elif isinstance(a, tuple):
        return tuple(map(wrap, a))
    elif isinstance(a, str):
        # this does not allocate memory to a pointer,
        # this only works if the destination to write to is an array and not a pointer
        return a.encode('utf-8')
    return a

def buildnew(f):
    return lambda *l: f(*wrap(l))

NULL = ffi.NULL

rcmfd2 = a()
for n in dir(rcmfd):
    setattr(rcmfd2, n, buildnew(getattr(rcmfd, n)))
rcmfd = rcmfd2
