"""
Define Material class and some material.

author : eablonet
date : 2019
"""


class Material(object):
    """Define a material."""
    def __init__(self):
        """Allocate values."""
        self.T = 0
        self.k = {}
        self.cp = {}
        self.rho = {}

    def __repr__(self):
        """Print action."""
        str = (
            "Values at {:d}°C\n\tcp : {:e}\n\t k : {:e}\n\t rho : {:e}" +
            "\n\t alpha : {:e}").format(
                self.T,
                self.get_cp(self.T),
                self.get_k(self.T),
                self.get_rho(self.T),
                self.get_alpha(self.T)
            )
        return str

    def get_cp(self, T=None):
        """Provide the specific capacity."""
        if T is not None:
            self.T = T
        Temps = sorted(self.cp)
        if self.T <= Temps[0]:
            return self.cp[Temps[0]]
        elif self.T >= Temps[-1]:
            return self.cp[Temps[-1]]
        else:
            i = 0
            while self.T > Temps[i]:
                i += 1
            p0 = (self.T - Temps[i-1]) / (Temps[i] - Temps[i-1])
            p1 = 1 - p0
            return p1*self.cp[Temps[i]] + p0*self.cp[Temps[i-1]]
        return

    def get_k(self, T=None):
        """Provide the conductivity."""
        if T is not None:
            self.T = T
        Temps = sorted(self.k)
        if self.T <= Temps[0]:
            return self.k[Temps[0]]
        elif self.T >= Temps[-1]:
            return self.k[Temps[-1]]
        else:
            i = 0
            while self.T > Temps[i]:
                i += 1
            p0 = (self.T - Temps[i-1]) / (Temps[i] - Temps[i-1])
            p1 = 1 - p0
            return p1*self.k[Temps[i]] + p0*self.k[Temps[i-1]]

    def get_rho(self, T=None):
        """Provide the density."""
        if T is not None:
            self.T = T
        Temps = sorted(self.rho)
        if self.T <= Temps[0]:
            return self.rho[Temps[0]]
        elif self.T >= Temps[-1]:
            return self.rho[Temps[-1]]
        else:
            i = 0
            while self.T > Temps[i]:
                i += 1
            p0 = (self.T - Temps[i-1]) / (Temps[i] - Temps[i-1])
            p1 = 1 - p0
            return p1*self.rho[Temps[i]] + p0*self.rho[Temps[i-1]]

    def get_alpha(self, T=None):
        """Calculate the value of diffusivity."""
        return self.get_k(T)/(self.get_rho(T)*self.get_cp(T))


class Water(Material):
    """Water properties."""

    def __init__(self):
        """Allocate values."""
        super().__init__()
        self.k = {
            0: 555.75e-3,
            10: 578.64e-3,
            20: 598.03e-3,
            30: 614.50e-3
        }  # source engineeringtoolbox.com
        self.rho = {
            0: 999.89,
            4: 999.95,
            10: 999.65,
            20: 998.19,
            30: 995.67
        }  # source engineeringtoolbox.com
        self.cp = {
            0: 4219.9,
            10: 4195.5,
            20: 4184.4,
            25: 4181.6,
            30: 4180.1
        }  # source engineeringtoolbox.com


class Ice(Material):
    """Define ice Material."""

    def __init__(self):
        """Allocate values."""
        super().__init__()
        self.k = {
            0: 2.22,
            -5: 2.25,
            -10: 2.30,
            -15: 2.34,
            -20: 2.39,
            -25: 2.45,
            -30: 2.50,
            -35: 2.57,
            -40: 2.63,
            -50: 2.76,
        }  # source engineeringtoolbox.com
        self.rho = {
            0: 916.2,
            -5: 917.5,
            -10: 918.9,
            -15: 919.4,
            -20: 919.4,
            -25: 919.6,
            -30: 920.0,
            -35: 920.4,
            -40: 920.8,
            -50: 921.6,
        }  # source engineeringtoolbox.com
        self.cp = {
            0: 2050,
            -5: 2027,
            -10: 2000,
            -15: 1972,
            -20: 1943,
            -25: 1913,
            -30: 1882,
            -35: 1851,
            -40: 1818,
            -50: 1751,
        }  # source engineeringtoolbox.com


class Copper(Material):
    """Define copper."""

    def __init__(self):
        """Allocate values."""
        super().__init__()
        self.rho = {
            0: 8940
        }
        self.k = {
            -73: 413,
            0: 401,
        }
        self.cp = {
            0: 390,
        }

class Aluminium(Material):
    def __init__(self):
        """Define values."""
        super().__init__()
        self.k = {
            -73: 237,
            0: 236,
            127: 240,
        }
        self.cp = {
            0: 910,
        }
        self.rho = {
            0: 2712,
        }





if __name__ == '__main__':
    co, wa, ic, al = Copper( ), Water(), Ice(), Aluminium()
    print(co)
    print(wa)
    print(ic)
    print(al)
