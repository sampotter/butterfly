import numpy as np

from matplotlib.colors import hsv_to_rgb

def complex_to_hsv(z, rmin=None, rmax=None, hue_start=90):
    amp = np.abs(z)
    if rmin is None:
        rmin = amp.min()
    if rmax is None:
        rmax = amp.max()

    # get amplidude of z and limit to [rmin, rmax]
    amp = np.where(amp < rmin, rmin, amp)
    amp = np.where(amp > rmax, rmax, amp)
    ph = np.angle(z, deg=1) + hue_start

    # HSV are values in range [0,1]
    h = (ph % 360)/360
    s = 0.85*np.ones_like(h)
    v = (amp - rmin)/(rmax - rmin)
    return hsv_to_rgb(np.dstack((h, s, v)))
