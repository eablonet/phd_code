# author : eablonet
# date # 2019

# load packages
# -------------
from numpy import linspace, zeros, arange
from numpy import sqrt, sum, exp, max
from numpy import count_nonzero, argmin
from scipy.special import erf, erfc
from numpy import pi
import warnings as warn
from matplotlib import pyplot as plt

from scipy import optimize as op

# locals packages
# ---------------
from packages.main import Material as ma
from packages.gui import ploter as pl
from packages.theory import schemes


class Stefan:
    """Stefan one dimensional analytical solution."""

    def __init__(self):
        """Init class."""
        # init solver data
        # ----------------
        self.solver = [0, 0, 1, 0]
        self.boundaries = [0, 0]
        self.boundValues = [-10, 0]
        self.__dilatation = False

        # init vectors
        # ------------
        self.mesh = Mesh()  # contain mesh, and time
        self.fields = FieldsVectors(self.mesh)  # contain alpha, T & zi

        # init data
        self.load_properties()
        self.set_time_condition()
        self.set_geometry()

    def load_properties(self):
        """Load physical parameters."""
        matlist = (
            ma.Water(), ma.Water('liquid'), ma.Water('solid'), ma.Copper(),
        )
        self.cp = zeros(4)
        self.k = zeros(4)
        self.rho = zeros(4)
        self.D = zeros(4)
        self.Tm = matlist[2].get_Tm()
        self.Lf = matlist[2].get_Lf()

        for i in range(4):
            mat = matlist[i]
            self.rho[i] = mat.get_rho()
            self.k[i] = mat.get_k()
            self.cp[i] = mat.get_cp()
            self.D[i] = mat.get_alpha()

    def set_boundaries(self, type=[0, 0], value=[-10, 0]):
        """Set boundaries condition, and value.

        input
        -----
            type: list
                Size of 2 int. First is the bottom, second is the top boundary.
                '0' for wall constant boundary, '1' for none flux boundary (
                dT/dz = 0), '2' for flux boundary (-kdT/dz = q).
            values : list
                Size of 2. Values to apply on bonddaries.
        """
        self.boundaries = type
        self.boundValues = value

    def set_geometry(self, H=[0, 1, 0], dz=[0, .1, 0], unit='mm'):
        """Define the problem geometry.

        input
        -----
            H: list
                SIze of 3. Init height of air, liquid and substrate.
            Nz: list
                Size of 3. Number of point in each domain.
            unit: str
                Height unit.
        """
        usi = {'mm': 1e-3, 'm': 1, 'cm': 1e-2}
        H = [h*usi[unit] for h in H]
        dz = [d*usi[unit] for d in dz]
        self.mesh.H = H
        self.mesh.zf = self.zf
        self.mesh.dz = dz

    def set_time_condition(self, ts=0, tend=10, dt=.1, auto=False):
        """Define time."""
        self.mesh.ts = ts
        self.mesh.tend = tend
        self.mesh.dt = dt
        self.autostop = auto

    @property
    def Ste(self):
        """Stefan number."""
        if self.solver[3] == 0:  # heat not solve in substrate
            return self.cp[2] * (self.Tm - self.boundValues[0]) / self.Lf

    @property
    def delta(self):
        """Return delta value."""
        return self.__delta

    @delta.setter
    def delta(self, val):
        self.__delta = val

    def _update_delta(self):
        if self.case == ('1:0:1:00'):  # ice-steady-wallwall
            delta = sqrt(self.Ste/2)
        elif self.case == ('1:0:1:20'):  # ice-steady-fluxwall
            delta = None
        elif self.case == ('1:0:1:22'):  # ice-steady-fluxflux
            None
        elif self.case == ('1:0:2:00'):  # ice-unsteady-wallwall
            def delta_method(delta):
                return delta*exp(delta**2)*erf(delta) - self.Ste/sqrt(pi)
            delta = op.fsolve(delta_method, .1)[0]
        elif self.case == ('1:0:2:20'):  # ice-unsteady-fluxwall
            None
        elif self.case == ('1:0:2:22'):  # ice-unsteady-fluxflux
            None
        elif self.case == '2:0:11:02':
            delta = None

        self.delta = delta

    @ property
    def tz0(self):
        """Return time to freeze initial height."""
        self._solver_check()
        self._case_selector()
        self._case_check()
        self._update_delta()

        # only ice cases
        # --------------
        if self.case == '1:0:1:00':  # ice-steady-wallwall
            return 1/(4*self.delta**2)*self.mesh.H[1]**2/self.D[2]
        elif self.case == '1:0:1:20':  # ice-steady-fluxwall
            return self.rho[2]*self.Lf / \
                (self.boundValues[0]*self.k[2])*self.mesh.H[1]
        elif self.case == '1:0:1:22':  # ice-steady-fluxflux
            None
        elif self.case == '1:0:2:00':  # ice-unsteady-wallwall
            return 1/(4*self.delta**2)*self.mesh.H[1]**2/self.D[2]
        elif self.case == '1:0:2:20':  # ice-unsteady-fluxwall
            None
        elif self.case == '1:0:2:22':  # ice-unsteady-fluxflux
            None

        # liquid-ice cases
        elif self.case == '2:0:11:00':
            None
        elif self.case == '2:0:11:01':
            None
        elif self.case == '2:0:11:02':
            if self.fields.zi is None:
                raise ValueError('Please to solve problem before calcul tz0')
            else:
                if max(self.fields.zi) < self.mesh.zf:
                    warn.warn('All the fields has not been solved.')
                return self.mesh.t[
                    argmin(abs(self.fields.zi - self.mesh.H[1]))
                ]
        elif self.case == '2:0:11:10':
            None
        elif self.case == '2:0:11:11':
            None
        elif self.case == '2:0:11:12':
            None
        elif self.case == '2:0:11:20':
            None
        elif self.case == '2:0:11:21':
            None
        elif self.case == '2:0:11:22':
            None

    @property
    def tf(self):
        """Return time to freezes all liquid.

        It's equal to tz0 in the absence of dilatation.
        """
        if self.dilatation is False:
            return self.tz0
        else:
            if self.boundValues == [0, 0]:
                return (self.rho[1]/self.rho[2])*self.tz0
            else:
                return (self.rho[1]/self.rho[2])**2*self.tz0

    @property
    def zf(self):
        """Return final height of ice."""
        if not self.dilatation:
            return self.mesh.H[1]
        else:
            return self.rho[1]/self.rho[2]*self.mesh.H[1]

    @property
    def dilatation(self):
        """Dilatation getter."""
        return self.__dilatation

    @dilatation.setter
    def dilatation(self, val):
        self.__dilatation = val
        self.mesh.zf = self.zf

    @property
    def solver(self):
        """Dilatation getter."""
        return self.__solver

    @solver.setter
    def solver(self, val):
        self.__solver = val
        self.n_phase = count_nonzero(self.solver)

    def _case_selector(self):
        """Choose the case of resolution."""
        self._solver_check()
        case = str(self.n_phase) + ':'
        situation = 1 if self.solver[3] != 0 else 0
        situation -= 1 if self.n_phase == 4 else 0
        case += str(situation) + ':'
        for s in self.solver:
            if s != 0:
                case += str(s)
        case += ':'
        for b in self.boundaries:
            case += str(b)
        self.case = case

    def _zi(self, n):
        """Calculate the solidification front position."""
        if self.case == ('1:0:1:00'):  # ice-steady-wallwall
            zi = 2*self.delta*sqrt(self.D[2] * self.mesh.t[n])
        elif self.case == ('1:0:1:20'):  # ice-steady-fluxwall
            zi = self.k[2]*self.boundValues[0] * self.mesh.t[n] / (self.rho[2] * self.Lf)
        elif self.case == ('1:0:1:22'):  # ice-steady-fluxflux
            None
        elif self.case == ('1:0:2:00'):  # ice-unsteady-wallwall
            zi = 2*self.delta*sqrt(self.D[2] * self.mesh.t[n])
        elif self.case == ('1:0:2:20'):  # ice-unsteady-fluxwall
            None
        elif self.case == ('1:0:2:22'):  # ice-unsteady-fluxflux
            None

        elif self.case == '2:0:11:02':
            if n == 0:  # __init__ zi to 0
                zi = 0
            elif n > 0:
                # load the previous time front position
                zn = self.fields.zi[n-1]
                if zn >= self.zf:
                    # on stop le calcul quand on atteint zf
                    zi = self.zf
                else:
                    # calculate new position
                    zi = sqrt(
                        zn**2 +
                        2*self.mesh.dt*(
                            self.D[2]*self.Ste +
                            self.boundValues[1]*zn/(self.rho[2]*self.Lf)
                        )
                    )
        return zi if zi <= self.zf else self.zf

    def _zl(self, n):
        if not self.dilatation:
            return self.mesh.H[1] - self.fields.zi[n]
        else:
            return self.mesh.H[1] - self.rho[2]/self.rho[1]*self.fields.zi[n]

    def _heat(self, j, n):
        """Return in a point at a time."""
        if self.case == ('1:0:1:00'):  # ice-steady-wallwall
            if self.mesh.z[j] == 0 or self.mesh.t[n] == 0:
                T = self.boundValues[0]
            elif self.mesh.z[j] <= self.fields.zi[n]:
                T = (
                    (self.Tm - self.boundValues[0]) *
                    self.mesh.z[j] / self.fields.zi[n] +
                    self.boundValues[0]
                )
            elif self.mesh.z[j] <= self.fields.zi[n] + self.fields.zl[n]:
                T = self.Tm
            else:
                return None
        elif self.case == ('1:0:1:20'):  # ice-steady-fluxwall
            if self.mesh.z[j] < self.fields.zi[n]:
                T = self.boundValues[0] * (self.mesh.z[j] - self.fields.zi[n]) + self.Tm
            elif self.mesh.z[j] <= self.fields.zi[n] + self.fields.zl[n]:
                T = self.Tm
            else:
                T = None
        elif self.case == ('1:0:1:22'):  # ice-steady-fluxflux
            None
        elif self.case == ('1:0:2:00'):  # ice-unsteady-wallwall
            if self.mesh.z[j] == 0 or self.mesh.t[n] == 0:
                T = self.Tm
            else:
                T = (
                    (
                        (self.Tm - self.boundValues[0]) * erf(
                            self.mesh.z[j] / (2*sqrt(self.D[2] * self.mesh.t[n]))
                            ) / erf(self.delta) + self.boundValues[0]
                    ) * (self.mesh.z[j] <= self.fields.zi[n]) +
                    self.Tm * (
                        (self.mesh.z[j] > self.fields.zi[n]) &
                        (self.mesh.z[j] <= self.fields.zi[n] + self.fields.zl[n])
                    )
                )

        elif self.case == ('1:0:2:20'):  # ice-unsteady-fluxwall
            None
        elif self.case == ('1:0:2:22'):  # ice-unsteady-fluxflux
            None

        elif self.case == '2:0:11:02':  # liquid-ice - steady - flux wall
            if self.mesh.z[j] == 0:
                T = self.Tm
            elif self.mesh.z[j] <= self.fields.zi[n]:
                T = (self.Tm - self.boundValues[0]) * self.mesh.z[j] / self.fields.zi[n] + self.boundValues[0]
            elif self.mesh.z[j] <= self.fields.zi[n] + self.fields.zl[n]:
                T = -self.boundValues[1] / self.k[1]*(self.mesh.z[j] - self.fields.zi[n]) + self.Tm
            else:
                T = None
        return T

    def _phase_function(self, z, t):
        None

    def _solver_check(self):
        if self.solver[2] == 0:  # ice unsolved
            raise ValueError(
                'Ice must be solved (1: steady, 2: transient)'
            )

        if self.n_phase == 2:
            if self.solver[0] != 0:  # try to solved air-ice
                raise ValueError(
                    "Only liquid-ice or ice-substrate for two phases could " +
                    "be solved"
                )

        if self.n_phase == 3:
            if self.solver[4] != 0 and self.solver[1] == 0:
                # with substrate & without liquid solved -> impossible
                raise ValueError(
                    'To solve air, liquid must be set up.' +
                    '(1: steady, 2: unsteady)'
                )

    def _case_check(self):
        if '1' in self.case.split(':')[3]:
            raise ValueError(
                'Non-FLux condition is not available in this case.'
            )

    def solve(self):
        """Solve the results."""
        self._case_selector()
        self._case_check()
        self._update_delta()

        if self.autostop:
            self.mesh.tend = self.tf

        self.fields.allocate_vectors()

        for n in range(self.mesh.Nt):
            self.fields.add_front_point(
                self._zi(n), n
            )
            self.fields.add_zl_point(
                self._zl(n), n
            )
            for j in range(self.mesh.Nz):
                self.fields.add_thermal_point(
                    self._heat(j, n),
                    j, n
                )


class Mesh:
    """Class containing mesh."""

    def __init__(self):
        """Init mashes Values."""
        # z-space
        self.__dz = None
        self.__H = None
        self.__zf = None  # final height of solidificaiton

        # t-space
        self.__dt = None
        self.__ts = None
        self.__tend = None

        # vectors allocation
        self.z = zeros(0)
        self.t = zeros(0)

    def _update_z(self):
        if (
            self.H is not None and
            self.dz is not None and
            self.zf is not None
        ):
            self.z = linspace(0, sum(self.H)-self.H[1]+self.zf, self.Nz)

    def _update_t(self):
        if (
            self.ts is not None and
            self.tend is not None and
            self.dt is not None
        ):
            self.t = linspace(self.ts, self.tend, self.Nt)

    @property
    def Nz(self):
        """Return the number of point for z discretization."""
        return int((sum(self.H)-self.H[1]+self.zf)/sum(self.dz))

    @property
    def Nt(self):
        """Return the number of point for t discretization."""
        return int((self.tend - self.ts)/self.dt)

    @property
    def dz(self):
        """Getter for dz."""
        return self.__dz

    @dz.setter
    def dz(self, val):
        self.__dz = val
        self._update_z()

    @property
    def H(self):
        """Getter for H."""
        return self.__H

    @H.setter
    def H(self, val):
        self.__H = val
        self._update_z()

    @property
    def zf(self):
        """Getter for H."""
        return self.__zf

    @zf.setter
    def zf(self, val):
        self.__zf = val
        self._update_z()

    @property
    def ts(self):
        """Getter for ts."""
        return self.__ts

    @ts.setter
    def ts(self, val):
        self.__ts = val
        self._update_t()

    @property
    def tend(self):
        """Getter for tend."""
        return self.__tend

    @tend.setter
    def tend(self, val):
        self.__tend = val
        self._update_t()

    @property
    def dt(self):
        """Getter for dt."""
        return self.__dt

    @dt.setter
    def dt(self, val):
        self.__dt = val
        self._update_t()


class FieldsVectors:
    """Class for Store and manipulate vectors."""

    def __init__(self, mesh):
        """Allocate vectors."""
        self.mesh = mesh

    def allocate_vectors(self):
        """Define vectors size."""
        self.zi = zeros(self.mesh.Nt)  # front location
        self.zl = zeros(self.mesh.Nt)  # front location
        self.alpha = zeros((self.mesh.Nz, self.mesh.Nt))  # phase function (0: liquid; 1: ice)
        self.T = zeros((self.mesh.Nz, self.mesh.Nt))

    def add_thermal_point(self, T, j, n):
        """Add point."""
        self.T[j, n] = T

    def add_front_point(self, zi, n):
        """Add a front point."""
        self.zi[n] = zi

    def add_zl_point(self, zl,  n):
        """Height of liquid versus time."""
        self.zl[n] = zl

    def add_alpha_point(self, alpha, j, n):
        """Add phase function point."""
        self.alpha[j, n] = alpha

    def plot(self):
        """Display ThermalField."""
        fig = plt.figure()
        cs = plt.contourf(
            self.mesh.t, self.mesh.z, self.T,
            cmap=plt.cm.bone, levels=10
        )
        plt.plot(self.mesh.t, self.zi, color='tab:red')
        plt.plot(self.mesh.t, self.zl+self.zi, color='tab:green')
        cbar = fig.colorbar(cs)
        cbar.ax.set_ylabel('Temperature', rotation=270)
        plt.show()

    def plot_liquid_air_temperature(self):
        """Display temperature at liquid-air interface."""
        j = argmin(abs(self.mesh.z - self.mesh.H[1]))
        fig = pl.Figure()
        ax = fig.add_ax(111)
        ax.plot(self.mesh.t, self.T[j, :])
        ax.xlabel('Time (s)')
        ax.ylabel('Temperature (°C)')
        fig.show()

if __name__ == '__main__':
    s = Stefan()

    s.solver = [0, 1, 1, 0]  # air, liq, ice, subs
    s.dilatation = False
    s.boundaries = [0, 2]  # bottom, top
    s.boundValues = [-10, -2000]  # bottom, top
    s.set_geometry(H=[0, 1.2, 0], dz=[0, .005, 0], unit='mm')  # air, liq, sub
    s.set_time_condition(ts=0, tend=20, dt=.01, auto=False)  # times

    s.solve()
    print('tf : ', s.tf)
    print('tz0 : ', s.tz0)
    s.fields.plot()
    s.fields.plot_liquid_air_temperature()
