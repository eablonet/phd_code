"""
List of method to generate quickly beatiful colors.

author: eablonet
date: #2019
"""
from matplotlib.pyplot import cm as cmap
from numpy import arange, array


def color_met(n):
    """Select beatiful colors in tab20c.

    input
    -----
        n: int
            number of color to generate. Maximum of 20 colors.
    """
    if n > 20:
        return IndexError('No more than 20pt')

    x = []
    a = n // 5
    b = n % 5
    for i in range(5):
        if i < b:
            x.extend(i*4 + arange(a+1))
        else:
            x.extend(i*4 + arange(a))
    c = cmap.tab20c(array(x).astype(int))
    return c

def color_grad(n, color='Blues'):
    if n > 20:
        raise IndexError('n must be less than 20')
    c = cmap.get_cmap(color, n)
    return c(arange(20))


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.figure()
    c = color_met(7)
    for i in range(7):
        plt.plot([0, 1], [i, i], c=c[i])

    plt.figure()
    c = color_grad(20)
    for i in range(20):
        plt.plot(0, i, marker='o', c=c[i])
    c = color_grad(20, 'Reds')
    for i in range(20):
        plt.plot(1, i, marker='o', c=c[-i])
    plt.show()
