"""
List of usefull methods.

author: eablonet
date: #2019

method
------
gradient(x, dt=1, order=2)
sliding_mean(x, size=3)

"""

from numpy import array, convolve, ones


def gradient(x, dt=1, order=2):
    """1D gradient.

    Compute the gradient dx/dt.

    inputs
    ------
        x : list / array
            Vector to dereivate
        dt : length beetween two points. Default 1.
        order : derivative order

    return
    ------
        dx/dt : derivate
    """
    if order == 2:
        mat = array([-1, 2, -1]) / 2

    x = convolve(x[::4], mat, 'same')

    return list(x) if type(x) is list else x
    # return a list if a list is given, else return an array


def sliding_mean(x, size=3):
    """Sliding average.

    Compute x(n-1)/3+x(n)/3+x(n+1)/3  # if size == 3

    inputs
    ------
        x : list / array
            Vector to dereivate
        dt : length beetween two points. Default 1.
        order : derivative order

    return
    ------
        dx/dt : derivate
    """
    x = convolve(x, ones((size,))/size, 'same')

    return list(x) if type(x) is list else x
