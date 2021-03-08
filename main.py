import iio
from rcmfd import rcmfd, NULL, STR, wrap
import numpy as np

im = iio.read('data/forged.tif').astype(np.float32)
ps = 8
tau = 1.
automatic = True
out = np.zeros(im.shape, dtype=np.float32)
im = np.ascontiguousarray(im.transpose(2,0,1))
c, h, w = im.shape
detection = rcmfd.perform_matching_py(c, im, w, h, ps, tau, automatic, out, False)
iio.write("test.png", out)

print("This image is a forgery: ", detection)
