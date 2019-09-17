# author: eablonet
# date : 2019

import matplotlib.pyplot as plt
import numpy as np


def table(T0, phi, plot=False):
    """Provide the Tint from tables."""
    Psat = [
        6.108, 6.112, 8.719, 10, 12.271, 15, 17.040, 20, 23.37, 25, 30,
        31.66, 40, 42.41
    ]
    Tsat = [
      0, 0.01, 5, 6.98, 10, 13.04, 15, 17.51, 20, 21.10, 24.10,
      25, 28.98, 30
    ]
    if phi > 1:
        phi /= 100

    c = 0
    while T0 > Tsat[c]:
        c += 1

    w = (T0 - Tsat[c]) / (Tsat[c+1] - Tsat[c])
    pv = phi * (w*Psat[c+1] + (1-w)*Psat[c])

    c = 0
    while pv > Psat[c]:
        c += 1

    w = (pv - Psat[c]) / (Psat[c+1] - Psat[c])
    Tint = (w*Tsat[c+1] + (1-w)*Tsat[c])

    if plot:
        plt.figure(figsize=[8, 4.5])
        plt.plot(Tsat, Psat, '--k')
        plt.plot(Tint, pv, 'or')
        plt.grid(True)
        plt.show()

    return Tint, pv


def clapeyron(T0, phi):
    """Provide Tint from Clapeyron formula.

    Asumption : Lv = Cst.
    """
    if phi > 1:
        phi /= 100

    r = 287

    Lv = 2453.6*1e3
    Tint = 1 / (1/T0 - np.log(phi)*r/Lv)

    return Tint


if __name__ == '__main__':
    Ta, phi = 20, 40
    Tint, pv = table(Ta, phi)
    print(Tint, pv)
    Tint_clap = clapeyron(Ta, phi)
    print(Tint_clap)

    To = []
    c = 0
    phi = range(10, 90, 1)
    Ta = [5, 10, 15, 20]
    for T in Ta:
        To.append([])
        for p in phi:
            Tint, pv = table(T, p)
            To[c].append(Tint)
        c += 1
    plt.figure(figsize=[8, 4.5])
    color = ['tab:blue', 'tab:red', 'tab:green', 'tab:orange']
    for i in range(4):
        plt.plot(
            phi, To[i], '-', color=color[i],
            label='Temperature ambiante : ' + str(Ta[i]) + '°C'
        )
    plt.grid(True)
    plt.xlabel(r'$\phi_a$ (%)', size=14)
    plt.ylabel(r'$T_{d}$', size=14)
    plt.ylim([0, 20])
    plt.legend(fancybox=True, shadow=True, fontsize=12)
    plt.title("Varaition de la température à l'interface liquide-air en fonction de l'humidité")
    plt.show()
