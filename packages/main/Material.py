"""
Define Material class and some material.

author : eablonet
date : 2019
"""


class Material(object):
    """Define a material."""

    def __init__(self, phase='solid'):
        """Allocate values."""
        self.T = 0  # temperature to select
        self. phase = phase  # phase to get

        self.k = {}  # thermal conductivity
        self.cp = {}  # specific heat
        self.rho = {}  # mass volumic

        self.Lf = 0  # fusion enthalpy
        self.Tm = 0  # melting temperature

        self.name = "None"

    def __repr__(self):
        """Print action."""
        str = (
                "{:s} ({:s})\n\t" +
                "Values at {:d}°C\n\t" +
                "cp : {:e}\n\t" +
                "k : {:e}\n\t" +
                "rho : {:e}\n\t" +
                "alpha : {:e}\n\t" +
                "Lf : {:e}\n\t" +
                "Tm : {:e}"
            ).format(
                self.name,
                self.phase,
                self.T,
                self.get_cp(self.T),
                self.get_k(self.T),
                self.get_rho(self.T),
                self.get_alpha(self.T),
                self.get_Lf(),
                self.get_Tm()
            )
        return str

    def set_phase(self, phase):
        """Change the phase to observe.

        input
        -----
            phase: string
                Choose the phase among:
                -liquid
                -solid
                -gaz
        """
        self.phase = phase

    def set_temperature(self, T):
        """Change the temperature of observation.

        input
        -----
            temp: float
                Temperature.
        """
        self.T = T

    def get_cp(self, T=None):
        """Provide the specific capacity."""
        if T is None:
            T = self.T

        Temps = sorted(self.cp[self.phase])
        if T <= Temps[0]:
            return self.cp[self.phase][Temps[0]]
        elif T >= Temps[-1]:
            return self.cp[self.phase][Temps[-1]]
        else:
            i = 0
            while T > Temps[i]:
                i += 1
            p0 = (T - Temps[i-1]) / (Temps[i] - Temps[i-1])

            return (1-p0)*self.cp[self.phase][Temps[i]] + \
                p0*self.cp[self.phase][Temps[i-1]]

    def get_k(self, T=None):
        """Provide the conductivity."""
        if T is None:
            T = self.T

        Temps = sorted(self.k[self.phase])
        if T <= Temps[0]:
            return self.k[self.phase][Temps[0]]
        elif T >= Temps[-1]:
            return self.k[self.phase][Temps[-1]]
        else:
            i = 0
            while T > Temps[i]:
                i += 1
            p0 = (T - Temps[i-1]) / (Temps[i] - Temps[i-1])

            return (1-p0)*self.k[self.phase][Temps[i]] + \
                p0*self.k[self.phase][Temps[i-1]]

    def get_rho(self, T=None):
        """Provide the density."""
        if T is None:
            T = self.T

        Temps = sorted(self.rho[self.phase])
        if T <= Temps[0]:
            return self.rho[self.phase][Temps[0]]
        elif T >= Temps[-1]:
            return self.rho[self.phase][Temps[-1]]
        else:
            i = 0
            while T > Temps[i]:
                i += 1
            p0 = (T - Temps[i-1]) / (Temps[i] - Temps[i-1])
            return (1-p0)*self.rho[self.phase][Temps[i]] + \
                p0*self.rho[self.phase][Temps[i-1]]

    def get_alpha(self, T=None):
        """Calculate the value of diffusivity."""
        if T is None:
            T = self.T

        return self.get_k(T)/(self.get_rho(T)*self.get_cp(T))

    def get_Tm(self, opt='celsius'):
        """Return the melting temperature."""
        usi = {'celsius': 0, 'kelvin': 273.16}
        return self.Tm + usi[opt]

    def get_Lf(self):
        """Return the enthalpy of fusion."""
        return self.Lf


class Water(Material):
    """Water properties."""

    def __init__(self, phase='solid'):
        """Allocate values."""
        super().__init__(phase)

        self.name = "Water"

        self.Lf = 335000
        self.Tm = 0

        self.k = {
            'liquid': {
                0: 555.75e-3,
                10: 578.64e-3,
                20: 598.03e-3,
                30: 614.50e-3
            },  # source engineeringtoolbox.com
            'solid': {
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
            },  # source engineeringtoolbox.com
        }
        self.rho = {
            'liquid': {
                0: 999.89,
                4: 999.95,
                10: 999.65,
                20: 998.19,
                30: 995.67
            },  # source engineeringtoolbox.com
            'solid': {
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
            },  # source engineeringtoolbox.com
        }
        self.cp = {
            'liquid': {
                0: 4219.9,
                10: 4195.5,
                20: 4184.4,
                25: 4181.6,
                30: 4180.1
            },  # source engineeringtoolbox.com
            'solid': {
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
            },  # source engineeringtoolbox.com
        }


class Copper(Material):
    """Define copper."""

    def __init__(self, phase='solid'):
        """Allocate values."""
        super().__init__(phase)

        self.name = "Copper"

        self.rho = {
            'solid': {
                0: 8940
            },  # source ??
        }
        self.k = {
            'solid': {
                -73: 413,
                0: 401,
            },  # source ??
        }
        self.cp = {
            'solid': {
                0: 390,
            },  # source ??
        }


class Aluminium(Material):
    """Define Aluminium."""

    def __init__(self, phase='solid'):
        """Define values."""
        super().__init__(phase)

        self.name = "Aluminium"

        self.k = {
            'solid': {
                -73: 237,
                0: 236,
                127: 240,
            },  # source ??
        }
        self.cp = {
            'solid': {
                0: 910,
            },  # source ??
        }
        self.rho = {
            'solid': {
                0: 2712,
            },  # source ??
        }


if __name__ == '__main__':
    co, wa, al = Copper(), Water(), Aluminium()
    print(wa)
    wa.set_phase('liquid')
    print(wa)
    print(co)
    print(al)
