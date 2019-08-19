"""
Stefan Solver Class.

----
The main class is the generic Stefan class.
Other class are depending on condition.
----
class
  Stefan : main

"""
import sys

from numpy import linspace, zeros, zeros_like, ones, arange, array
from numpy import sqrt, arccos, exp
from scipy.special import erf, erfc
from numpy import pi

from matplotlib import pyplot as plt

from scipy import optimize as op

try:
    from packages.main import Material as ma
except ModuleNotFoundError:
    import Material as ma


class Stefan():
    """Help to select the good Stefan solution."""

    def __init__(self):
        """Init the class."""
        self.dimension = 1
        self.set_ice_model()  # init ice model
        self.set_liquid_model()  # init liquid model
        self.set_air_model()  # init air model
        self.set_wall_model()  # init wall model
        self.set_dilatation()  # ini dilatation

    def set_ice_model(self, mode=1):
        """Change ice model.

        inputs
        ------
            mode : int
                1, 2. Default 1.
                1 : stationary diffusion
                2 : instationary diffusion
        """
        self.ice_model = mode

    def set_liquid_model(self, mode=0):
        """Change ice model.

        inputs
        ------
            mode : int
                0, 1, 2, 3. Default 0.
                0 : unsolved
                1 : stationary diffusion
                2 : instationary diffusion
                3 : instationary advection-diffusion
        """
        self.liquid_model = mode

    def set_air_model(self, mode=0):
        """Change ice model.

        inputs
        ------
            mode : int
                0, 1, 2. Default 0.
                0 : unsolved
                1 : stationary diffusion
                2 : instationary diffusion
        """
        self.air_model = mode

    def set_dilatation(self, val=False):
        """Set the dilatation solution.

        Caution, the dilatation is 1D. So:
            Final z:
            zf = rho_l/rho_i z0 else zf = z0
            Liquid height time evolution:
            zl(t) = z0 - rho_i/rho_l zi(t) else zl(t) = z0 - zi(t)
            If advection is solved, the dilatation modify the velocity in the
            liquid. Note that is there is no dilatation advection is null in
            the liquid.

        inputs
        ------
            val : bool
                Default False. True to solve take care of dilatation.
        """
        self.dilatation = val

    def get_theorical(self):
        """Return goog theoritical model."""
        case = {
            1000: 'Ste_1d_1c_1p_qs',
            2000: 'Ste_1d_1c_1p_is',
            1200: 'Ste_1d_1c_2p_qsis',
            2200: 'Ste_1d_1c_2p_isis',
        }
        cons = (
            str(self.ice_model) +
            str(self.liquid) +
            str(self.air_model) +
            str(0)
        )
        return case[cons](self.dilatation)


class Ste_1d_1c_1p_qs(object):
    """Solve Stefan quasi-static, only in solid phase.

    1c : one component : water
    1p : one phase : ice

    Heat transfer with phase chage solution of :
        1. ∆T = 0
        2. when z->0: T = Tw
        3. when z-> z_i:
            T(z_i(t)) = Tm
            rho_i Lf dz/dt = k dT/dz

    Parameters
    ----------
        rho: float
            material density
        k: float
            thermal conductivity
        cp: float
            specific capacity
        alpha: float
            diffusivity
        Lf: float
            Fusion enthalpy
        Lf: float
            Fusion enthalpy

    Methods
    -------
        import_material

    """

    def __init__(self, dilatation=False):
        """Allocate variables.

        input
        -----
            mat: string
                Material to solidify. Default water.

        """
        # init physical paramaters
        # ------------------------
        self.rho = {'solid': None}
        self.k = None
        self.cp = None
        self.alpha = None
        self.Lf = None

        # init dilatation parameters
        # --------------------------
        self.dilatation = dilatation
        if dilatation:
            self.rho['liquid'] = None

        # init thermal paramaters
        # -----------------------
        self.Tw = None  # wall temperature
        self.Tm = None  # melting temperature
        self.T = None  # temperature vector in ice

        # init time
        # ---------
        self.time = None

        # init geometry values
        # --------------------
        self.z0 = None  # initiale height of liquid
        self.zf = None  # final time of solidification
        self.z = None  # space vector

        # init output values
        # ------------------
        self.zi = None  # ice height
        self.zl = None  # liquid height  # ndrl z_tot = zl+zi
        self.Ste = None  # stefan number
        self.delta = None  # delta value
        self.tz0 = None  # time to reach z0
        self.tf = None  # time to freeze all the liquid, ...
        # equal to tz0 without dilatation

    def __repr__(self):
        """Represent function."""
        return (
            'One dimensional steady monophasic solution of Stefan Problem.'
        )

    def set_dilatation(self, val):
        """Enable/disable dilatation.

        inputs
        ------
            val : bool
                True or False. Default is False.
        """
        self.dilatation = val

    def set_material(self, material='Water'):
        """Import the material properties."""
        mat = {'Water': ma.Water(), }
        mat = mat[material]
        mat.set_temperature(mat.Tm)

        self.rho['solid'] = mat.get_rho()
        self.k = mat.get_k()
        self.cp = mat.get_cp()
        self.alpha = mat.get_alpha()

        self.Lf = mat.get_Lf()
        self.Tm = mat.get_Tm()

        if self.dilatation:
            mat.set_phase('liquid')
            self.rho['liquid'] = mat.get_rho()

    def set_geometry(self, z0, zmax=None, Nz=100):
        """Import geometry info, create geometry vector.

        inputs
        ------
            z0 : float
                initial height of liquid
            zmax : float
                height in which calculate thermal solution. For z > z0,
                it's assume that it's the same material than initial, but there
                is not phase transition.

        """
        self.z0 = z0
        if self.dilatation:
            if self.rho['liquid'] is not None:
                self.zf = self.rho['liquid']/self.rho['solid']*z0
            else:
                raise ValueError('Please to load material before geometry')
        else:
            self.zf = z0

        if zmax is None:
            if self.dilatation:
                zmax = self.zf
            else:
                zmax = z0
        self.z = linspace(0, zmax, Nz)

    def set_thermic(self, T):
        """Set the wall temperature.

        inputs
        ------
            T : float
                Is the wall temperature
        """
        self.Tw = T

    def set_time(self, tmax=1, Nt=100, auto=False):
        """Import geometry info, create geometry vector.

        inputs
        ------
            tf: float
                total time to make the computation
            Nt: int
                number of point fro discretization. Default 100.
            auto: bool
                Auto calculate the total time of solidification. Default False.
                If True, it need z0 and thermics. So please run set_geometry
                and set_thermic before set_time.
                In this case tf is ignored.
        """
        if auto:
            if self.z0 is None or self.Tw is None:
                raise ValueError(
                    'Please run set_geometry and set_thermic to use auto ' +
                    'time calculation'
                )
            else:
                if self.dilatation:
                    zf = self.zf
                else:
                    zf = self.z0
                tmax = zf**2*self.Lf / (
                    2*self.alpha*self.cp*(self.Tm - self.Tw)
                )
        self.time = linspace(0, tmax, Nt)

    def solve(self):
        """Solve front position, and time field.

        Also calculate: delta, Stefan and tz0.
        """
        # init temperature field with Tm
        # ------------------------------
        self.T = self.Tm*ones([len(self.z), len(self.time)])

        # zi solution
        # -----------
        zi = sqrt(
            2*self.k*(self.Tm - self.Tw)*self.time /
            (self.rho['solid'] * self.Lf)
        )
        if self.dilatation:
            zf = self.zf
        else:
            zf = self.z0
        self.zi = zi*(zi < zf) + zf*(zi >= zf)

        # zl calcul
        # ---------
        # mass conservation : \rho_i == \rho_l
        if self.dilatation:
            self.zl = self.z0 - self.rho['solid']/self.rho['liquid']*self.zi
        else:
            self.zl = self.z0 - self.zi

        # calc tz0, tf, Stefan and delta
        # ------------------------------
        self.tz0 = self.z0**2*self.Lf / (
            2*self.alpha*self.cp*(self.Tm - self.Tw)
        )
        self.tf = self.zf**2*self.Lf / (
            2*self.alpha*self.cp*(self.Tm - self.Tw)
        )
        self.Ste = self.cp * (self.Tm - self.Tw) / self.Lf
        self.delta = sqrt(self.Ste/2)

        # calculate T for each time
        # -------------------------
        # si z < zi, then we are in the solid (diffusion), else we are in the
        # liquid (T=cste).
        for t in range(1, len(self.time)):
            self.T[:, t] = (
                (self.Tm - self.Tw)*self.z/self.zi[t] + self.Tw
            ) * (self.z < self.zi[t])

    def get_Ste(self):
        """Return Stefan value."""
        return self.Ste

    def get_delta(self):
        """Return delta value."""
        return self.delta

    def get_tz0(self):
        """Return time to freeze z0 value."""
        return self.tz0

    def get_tf(self):
        """Return time to frezze all the liquid value.

        This is equal to tz0 is there is no dilatation.
        """
        return self.tf

    def get_zf(self):
        """Return final time of solidification."""
        return self.zf

    def get_front(self, adim=False):
        """Return the solid position over time.

        inputs
        ------
            adim : bool
                if True, return the adimentional temperature,
                if False, return the dimensional one.

        """
        return self.zi if not adim else self.zi/self.z0

    def get_tempfield(self, adim=False):
        """Return the temperature field.

        return = T(z, t)

        Parameters
        ----------
        adim : bool
            if True, return the adimentional temperature,
            if False, return the dimensional one.

        """
        if adim:
            T = (self.T - self.Tw)/(self.Tm - self.Tw)
        else:
            T = self.T
        return T

    def get_time(self):
        """Return time vector."""
        return self.time

    def get_mesh(self):
        """Return space vector."""
        return self.z

    def get_liquid_height(self):
        """Return the liquid height of the drop over time."""
        return self.zl

    def get_total_height(self):
        """Return the total height of the drop over time."""
        return self.zi + self.zl

    def plot_front(self, show=True, color=None):
        """Graphical solution of stefan problem."""
        plt.figure(figsize=[8, 4.5])
        if color is None:
            color = 'tab:blue'

        plt.plot(self.time, self.zi*1e3, '-', color=color)

        plt.plot([self.tf, self.tf], [0, self.zf*1e3], '.--', c=color)
        plt.plot([0, self.tf], [self.zf*1e3, self.zf*1e3], '.--', c=color)

        if self.dilatation:
            plt.plot(
                [self.tz0, self.tz0], [0, self.z0*1e3],
                '.--', c=color
            )
            plt.plot(
                [0, self.tz0], [self.z0*1e3, self.z0*1e3],
                '.--', c=color
            )

        plt.xlabel('time (s)')
        plt.ylabel('z (mm)')
        plt.grid(True)

        if show:
            plt.show()

    def plot_heat_field(self):
        """Graphical representation of the thermal field."""
        fig = plt.figure(figsize=[8, 4.5])
        ax = plt.subplot(1, 1, 1)
        cs = ax.contourf(
            self.time, self.z*1e3, self.T,
            cmap=plt.cm.bone, levels=10
        )
        ax.plot(self.time, (self.zi + self.zl)*1e3, '--', c='k')
        ax.plot(self.time, self.zi*1e3, '--', c='tab:blue')
        cbar = fig.colorbar(cs)
        cbar.ax.set_ylabel('Temperature', rotation=270)


        ax.set_xlabel('temps (s)')
        ax.set_ylabel('z (mm)')
        plt.show()

class Ste_1d_1c_1p_is(Ste_1d_1c_1p_qs):
    """Solve Stefan instationny, only in solid phase.

    1c : one component : Water
    1p : one phase : solid

    Heat transfer with phase chage solution of :
        1. dT/dt = D ∆T
        2. when z->0: T = Tw
        3. when z-> z_i:
            T(z_i(t)) = Tm
            rho_i Lf dz/dt = k dT/dz
    """

    def __init__(self, dilatation=False):
        """Allocate.

        Nothing new from Stefan_1c_1p_qs. Only the solve method changes.
        """
        super().__init__(dilatation)

    def set_time(self, tmax=1, Nt=100, auto=False):
        """Import geometry info, create geometry vector.

        inputs
        ------
            tf: float
                total time to make the computation
            Nt: int
                number of point fro discretization. Default 100.
            auto: bool
                Auto calculate the total time of solidification. Default False.
                If True, it need z0 and thermics. So please run set_geometry
                and set_thermic before set_time.
                In this case tf is ignored.
        """
        if auto:
            if self.z0 is None or self.Tw is None:
                raise ValueError(
                    'Please run set_geometry and set_thermic to use auto ' +
                    'time calculation'
                )
            else:
                def delta_method(delta):
                    Ste = self.cp * (self.Tm - self.Tw) / self.Lf
                    return delta*exp(delta**2)*erf(delta) - Ste/sqrt(pi)
                delta = op.fsolve(delta_method, .1)[0]
                # .1 because we expect a small delta
                if self.dilatation:
                    zf = self.zf
                else:
                    zf = self.z0
                tmax = zf**2 / (4*delta**2*self.alpha)

        self.time = linspace(0, tmax, Nt)

    def solve(self):
        """Solve front position, and time field.

        Also calculate: delta, Stefan and tz0.
        """
        # enable/disable dilatation
        if self.dilatation:
            zf = self.zf
        else:
            zf = self.z0

        # init temperature field with Tm
        # ------------------------------
        self.T = self.Tm*ones([len(self.z), len(self.time)])

        # calc tz0, Stefan and delta
        # --------------------------
        self.Ste = self.cp * (self.Tm - self.Tw) / self.Lf

        # calculate delta
        # ---------------
        def delta_method(delta):
            return delta*exp(delta**2)*erf(delta) - self.Ste/sqrt(pi)
        self.delta = op.fsolve(delta_method, .1)[0]
        # .1 because we expect a small delta

        # calc tz0
        # --------
        self.tz0 = self.z0**2 / (4*self.delta**2*self.alpha)

        # calc tf
        # --------
        self.tf = self.zf**2 / (4*self.delta**2*self.alpha)

        # zi solution
        # -----------
        zi = 2 * self.delta*sqrt(self.alpha*self.time)
        self.zi = zi*(zi < zf) + zf*(zi >= zf)

        # zl calcul
        # ---------
        # mass conservation : \rho_i == \rho_l
        if self.dilatation:
            self.zl = self.z0 - self.rho['solid']/self.rho['liquid']*self.zi
        else:
            self.zl = self.z0 - self.zi

        # calculate T for each time
        # -------------------------
        # si z < zi, then we are in the solid (diffusion), else we are in the
        # liquid (T=cste).
        for t in range(1, len(self.time)):
            self.T[:, t] = (
                (self.Tm - self.Tw) * erf(
                    self.z/(2*sqrt(self.alpha*self.time[t]))
                ) / erf(self.delta) + self.Tw
            ) * (self.z < self.zi[t])


class Ste_1d_1c_2p_qsis(Ste_1d_1c_1p_qs):
    """Solve Stefan instationny, only in solid phase.

    1c : one component : Water
    2p : two phase : solid - liquid
    qsis : solid quasi static, liquid instationary

    Heat transfer with phase chage solution of :
        Liquid :
        1. dT/dt = D ∆T
        2. when z->: T = T_inf
        3. when z-> z_i:
            T(z_i(t)) = Tm
            rho_i Lf dz/dt = q_l - q_s
        Solid :
        1. dT/dt = 0
        2. when z->0: T = Tw
        3. when z-> z_i:
            T(z_i(t)) = Tm
            rho_i Lf dz/dt = q_l - q_s
    """

    def __init__(self, dilatation=False):
        """Allocate variables.

        input
        -----
            mat: string
                Material to solidify. Default water.

        """
        super().__init__(dilatation)
        # init physical paramaters
        # ------------------------
        self.rho = {'liquid': None, 'solid': None}
        self.k = {'liquid': None, 'solid': None}
        self.cp = {'liquid': None, 'solid': None}
        self.alpha = {'liquid': None, 'solid': None}

        # init thermal paramaters
        # -----------------------
        self.Tinf = None  # temperature in liquid far from wall (== t=0)

        # output data
        # -----------
        self.Ste = {'solid': None, 'liquid': None}

    def __repr__(self):
        """Represent function."""
        return (
            'One dimensional steady monophasic solution of Stefan Problem.'
        )

    def set_material(self, material='Water'):
        """Import the material properties."""
        mat = {'Water': ma.Water(), }
        mat = mat[material]
        mat.set_temperature(mat.Tm)
        self.Lf = mat.get_Lf()
        self.Tm = mat.get_Tm()

        for k in ['liquid', 'solid']:
            mat.set_phase(k)
            self.rho[k] = mat.get_rho()
            self.k[k] = mat.get_k()
            self.cp[k] = mat.get_cp()
            self.alpha[k] = mat.get_alpha()

    def set_thermic(self, T=[-15, 10]):
        """Set the wall temperature.

        inputs
        ------
            T : vector, list, array... of size 2
                T[0] is the wall temperature, and T[1] is the infinite one.
        """
        self.Tw = T[0]
        self.Tinf = T[1]

    def set_time(self, tmax=30, Nt=100, auto=False):
        """Import geometry info, create geometry vector.

        inputs
        ------
            tf: float
                total time to make the computation
            Nt: int
                number of point fro discretization. Default 100.
            auto: bool
                Auto calculate the total time of solidification. Default False.
                If True, it need z0 and thermics. So please run set_geometry
                and set_thermic before set_time.
                In this case tf is ignored.
        """
        if auto:
            if self.z0 is None or self.Tw is None:
                raise ValueError(
                    'Please run set_geometry and set_thermic to use auto ' +
                    'time calculation'
                )
            else:
                # calculate delta
                # ---------------
                def delta_method(delta):
                    Stes = self.cp['solid'] * \
                        (self.Tm - self.Tw) / self.Lf
                    Stel = self.cp['liquid'] * \
                        (self.Tinf - self.Tm) / self.Lf
                    f = delta**2 - \
                        Stes/2 + \
                        delta*self.rho['liquid'] / self.rho['solid'] * \
                        sqrt(self.alpha['liquid'] / self.alpha['solid']) * \
                        Stel/sqrt(pi) * \
                        exp(
                            -delta**2 *
                            self.alpha['solid']/self.alpha['liquid']
                        ) / \
                        erfc(
                            delta * sqrt(
                                self.alpha['solid']/self.alpha['liquid']
                            )
                        )
                    return f
                delta = op.fsolve(delta_method, .1)[0]
                # .1 because we expect a small delta
                if self.dilatation:
                    zf = self.zf
                else:
                    zf = self.z0
                tmax = zf**2 / (4*delta**2*self.alpha['solid'])
                # .1 because we expect a small delta

        self.time = linspace(0, tmax, Nt)

    def solve(self):
        """Solve front position, and time field.

        Also calculate: delta, Stefan and tz0.
        """
        # enable / disable dilatation
        # ---------------------------
        if self.dilatation:
            zf = self.zf
        else:
            zf = self.z0

        # init temperature field with Tm
        # ------------------------------
        self.T = zeros([len(self.z), len(self.time)])

        # calc tz0, Stefan and delta
        # --------------------------
        self.Ste['solid'] = self.cp['solid'] * \
            (self.Tm - self.Tw) / self.Lf
        self.Ste['liquid'] = self.cp['liquid'] * \
            (self.Tinf - self.Tm) / self.Lf

        # calculate delta
        # ---------------
        def delta_method(delta):
            f = delta**2*erfc(
                    delta*sqrt(self.alpha['solid']/self.alpha['liquid'])
                ) - \
                self.Ste['solid']/2 * erfc(
                    delta*sqrt(self.alpha['solid']/self.alpha['liquid'])
                ) + \
                delta*self.rho['liquid'] / self.rho['solid'] * \
                sqrt(self.alpha['liquid'] / self.alpha['solid']) * \
                self.Ste['liquid'] / sqrt(pi) * \
                exp(
                    -delta**2*self.alpha['solid']/self.alpha['liquid']
                )

            return f
        self.delta = op.fsolve(delta_method, .1)[0]
        # .1 because we expect a small delta

        # calc tz0
        # --------
        self.tz0 = self.z0**2 / (4*self.delta**2*self.alpha['solid'])

        # calc tf
        # --------
        self.tf = self.zf**2 / (4*self.delta**2*self.alpha['solid'])

        # zi solution
        # -----------
        zi = 2 * self.delta*sqrt(self.alpha['solid']*self.time)
        self.zi = zi*(zi < zf) + zf*(zi >= zf)

        # zl calcul
        # ---------
        # mass conservation : \rho_i == \rho_l
        if self.dilatation:
            self.zl = self.z0 - self.rho['solid']/self.rho['liquid']*self.zi
        else:
            self.zl = self.z0 - self.zi

        # calculate T for each time
        # -------------------------
        # si z < zi, then we are in the solid (diffusion), else we are in the
        # liquid (T=cste).
        for t in range(1, len(self.time)):
            self.T[self.z < self.zi[t], t] = (
                (self.Tm - self.Tw) *
                self.z[self.z < self.zi[t]] / self.zi[t] +
                self.Tw
            )
            self.T[self.z > self.zi[t], t] = (
                (self.Tm - self.Tinf) * erfc(
                    self.z[self.z > self.zi[t]] /
                    (2*sqrt(self.alpha['liquid']*self.time[t]))
                ) /
                erfc(
                    self.delta*sqrt(self.alpha['solid']/self.alpha['liquid'])
                ) +
                self.Tinf
            )
            self.T[self.z == self.zi[t], t] = self.Tm


class Ste_1d_1c_2p_isis(Ste_1d_1c_2p_qsis):
    """Solve Stefan instationnary for both liquid and solid.

    1c : one component : Water
    2p : two phase : solid - liquid
    isis : solid and liquid both instationnary.

    Heat transfer with phase chage solution of :
        Liquid :
        1. dT/dt = D ∆T
        2. when z->: T = T_inf
        3. when z-> z_i:
            T(z_i(t)) = Tm
            rho_i Lf dz/dt = q_l - q_s
        Solid :
        1. dT/dt = D ∆T
        2. when z->0: T = Tw
        3. when z-> z_i:
            T(z_i(t)) = Tm
            rho_i Lf dz/dt = q_l - q_s
    """

    def set_time(self, tmax=30, Nt=100, auto=False):
        """Import geometry info, create geometry vector.

        inputs
        ------
            tf: float
                total time to make the computation
            Nt: int
                number of point fro discretization. Default 100.
            auto: bool
                Auto calculate the total time of solidification. Default False.
                If True, it need z0 and thermics. So please run set_geometry
                and set_thermic before set_time.
                In this case tf is ignored.
        """
        if auto:
            if self.z0 is None or self.Tw is None:
                raise ValueError(
                    'Please run set_geometry and set_thermic to use auto ' +
                    'time calculation'
                )
            else:
                # calculate delta
                # ---------------
                def delta_method(delta):
                    Stes = self.cp['solid'] * \
                        (self.Tm - self.Tw) / self.Lf
                    Stel = self.cp['liquid'] * \
                        (self.Tinf - self.Tm) / self.Lf
                    f = delta*erf(delta)*exp(delta**2) - \
                        Stes/sqrt(pi) + \
                        self.rho['liquid'] / self.rho['solid'] * \
                        sqrt(self.alpha['liquid'] / self.alpha['solid']) * \
                        Stel/sqrt(pi) * \
                        exp(
                            delta**2*(
                                1-self.alpha['solid']/self.alpha['liquid']
                            )
                        ) / \
                        erfc(
                            delta*sqrt(
                                self.alpha['solid']/self.alpha['liquid']
                            )
                        ) * \
                        erf(delta)
                    return f
                delta = op.fsolve(delta_method, .1)[0]
                # .1 because we expect a small delta

                if self.dilatation:
                    zf = self.zf
                else:
                    zf = self.z0
                tmax = zf**2 / (4*delta**2*self.alpha['solid'])
                # .1 because we expect a small delta

        self.time = linspace(0, tmax, Nt)

    def solve(self):
        """Solve front position, and time field.

        Also calculate: delta, Stefan and tz0.
        """
        # enable / disable dilatation
        # ---------------------------
        if self.dilatation:
            zf = self.zf
        else:
            zf = self.z0

        # init temperature field with Tm
        # ------------------------------
        self.T = zeros([len(self.z), len(self.time)])

        # calc tz0, Stefan and delta
        # --------------------------
        self.Ste['solid'] = self.cp['solid'] * \
            (self.Tm - self.Tw) / self.Lf
        self.Ste['liquid'] = self.cp['liquid'] * \
            (self.Tinf - self.Tm) / self.Lf

        # calculate delta
        # ---------------
        def delta_method(delta):
            f = delta*erf(delta)*exp(delta**2) - \
                self.Ste['solid']/sqrt(pi) + \
                self.rho['liquid'] / self.rho['solid'] * \
                sqrt(self.alpha['liquid'] / self.alpha['solid']) * \
                self.Ste['liquid']/sqrt(pi) * \
                exp(delta**2*(1-self.alpha['solid']/self.alpha['liquid'])) / \
                erfc(delta*sqrt(self.alpha['solid']/self.alpha['liquid'])) * \
                erf(delta)
            return f
        self.delta = op.fsolve(delta_method, .1)[0]
        # .1 because we expect a small delta

        # calc tz0
        # --------
        self.tz0 = self.z0**2 / (4*self.delta**2*self.alpha['solid'])

        # calc tf
        # --------
        self.tf = self.zf**2 / (4*self.delta**2*self.alpha['solid'])

        # zi solution
        # -----------
        zi = 2 * self.delta*sqrt(self.alpha['solid']*self.time)
        self.zi = zi*(zi < zf) + zf*(zi >= zf)

        # zl calcul
        # ---------
        # mass conservation : \rho_i == \rho_l
        if self.dilatation:
            self.zl = self.z0 - self.rho['solid']/self.rho['liquid']*self.zi
        else:
            self.zl = self.z0 - self.zi

        # calculate T for each time
        # -------------------------
        # si z < zi, then we are in the solid (diffusion), else we are in the
        # liquid (T=cste).
        for t in range(1, len(self.time)):
            self.T[self.z < self.zi[t], t] = (
                (self.Tm - self.Tw) *
                erf(
                    self.z[self.z < self.zi[t]] / (
                        2*sqrt(self.alpha['solid']*self.time[t])
                    )
                ) /
                erf(self.delta) +
                self.Tw
            )
            self.T[self.z > self.zi[t], t] = (
                (self.Tm - self.Tinf) * erfc(
                    self.z[self.z > self.zi[t]] /
                    (2*sqrt(self.alpha['liquid']*self.time[t]))
                ) /
                erfc(
                    self.delta*sqrt(self.alpha['solid']/self.alpha['liquid'])
                ) +
                self.Tinf
            )
            self.T[self.z == self.zi[t], t] = self.Tm


class Ste_1d_1c_2p_isad(Ste_1d_1c_2p_isis):
    """Solve Stefan instationnary for both liquid and solid.

    1c : one component : Water
    2p : two phase : solid - liquid
    isad : solid instationnary and advection in liquid.

    Heat transfer with phase chage solution of :
        Liquid :
        1. dT/dt + u grad(T) = D ∆T
        2. when z->: T = T_inf
        3. when z-> z_i:
            T(z_i(t)) = Tm
            rho_i Lf dz/dt = q_l - q_s
        Solid :
        1. dT/dt = D ∆T
        2. when z->0: T = Tw
        3. when z-> z_i:
            T(z_i(t)) = Tm
            rho_i Lf dz/dt = q_l - q_s
    """

    def set_time(self, tmax=30, Nt=100, auto=False):
        """Import geometry info, create geometry vector.

        inputs
        ------
            tf: float
                total time to make the computation
            Nt: int
                number of point fro discretization. Default 100.
            auto: bool
                Auto calculate the total time of solidification. Default False.
                If True, it need z0 and thermics. So please run set_geometry
                and set_thermic before set_time.
                In this case tf is ignored.
        """
        if auto:
            if self.z0 is None or self.Tw is None:
                raise ValueError(
                    'Please run set_geometry and set_thermic to use auto ' +
                    'time calculation'
                )
            else:
                # calculate delta
                # ---------------
                def delta_method(delta):
                    Stes = self.cp['solid'] * \
                        (self.Tm - self.Tw) / self.Lf
                    Stel = self.cp['liquid'] * \
                        (self.Tinf - self.Tm) / self.Lf
                    f = delta*erf(delta)*exp(delta**2) - \
                        Stes/sqrt(pi) + \
                        self.rho['liquid'] / self.rho['solid'] * \
                        sqrt(self.alpha['liquid'] / self.alpha['solid']) * \
                        Stel/sqrt(pi) * \
                        exp(
                            delta**2*(
                                1 -
                                self.rho['solid']**2/self.rho['liquide']**2 *
                                self.alpha['solid']/self.alpha['liquid']
                            )
                        ) / \
                        erfc(
                            delta *
                            self.rho['solid']/self.rho['liquide'] *
                            sqrt(
                                self.alpha['solid']/self.alpha['liquid']
                            )
                        ) * \
                        erf(delta)
                    return f
                delta = op.fsolve(delta_method, .1)[0]
                # .1 because we expect a small delta

                if self.dilatation:
                    zf = self.zf
                else:
                    zf = self.z0
                tmax = zf**2 / (4*delta**2*self.alpha['solid'])
                # .1 because we expect a small delta

        self.time = linspace(0, tmax, Nt)

    def solve(self):
        """Solve front position, and time field.

        Also calculate: delta, Stefan and tz0.
        """
        # enable / disable dilatation
        # ---------------------------
        if self.dilatation:
            zf = self.zf
        else:
            zf = self.z0

        # init temperature field with Tm
        # ------------------------------
        self.T = zeros([len(self.z), len(self.time)])

        # calc tz0, Stefan and delta
        # --------------------------
        self.Ste['solid'] = self.cp['solid'] * \
            (self.Tm - self.Tw) / self.Lf
        self.Ste['liquid'] = self.cp['liquid'] * \
            (self.Tinf - self.Tm) / self.Lf

        # calculate delta
        # ---------------
        def delta_method(delta):
            f = delta*erf(delta)*exp(delta**2) - \
                self.Ste['solid']/sqrt(pi) + \
                self.rho['liquid'] / self.rho['solid'] * \
                sqrt(self.alpha['liquid'] / self.alpha['solid']) * \
                self.Ste['liquid']/sqrt(pi) * \
                exp(
                    delta**2*(
                        1 -
                        self.rho['solid']**2/self.rho['liquid']**2 *
                        self.alpha['solid']/self.alpha['liquid']
                    )
                ) / \
                erfc(
                    delta *
                    self.rho['solid']/self.rho['liquid'] *
                    sqrt(self.alpha['solid']/self.alpha['liquid'])
                ) * \
                erf(delta)
            return f
        self.delta = op.fsolve(delta_method, .1)[0]
        # .1 because we expect a small delta

        # calc tz0
        # --------
        self.tz0 = self.z0**2 / (4*self.delta**2*self.alpha['solid'])

        # calc tf
        # --------
        self.tf = self.zf**2 / (4*self.delta**2*self.alpha['solid'])

        # zi solution
        # -----------
        zi = 2 * self.delta*sqrt(self.alpha['solid']*self.time)
        self.zi = zi*(zi < zf) + zf*(zi >= zf)

        # zl calcul
        # ---------
        # mass conservation : \rho_i == \rho_l
        if self.dilatation:
            self.zl = self.z0 - self.rho['solid']/self.rho['liquid']*self.zi
        else:
            self.zl = self.z0 - self.zi

        # calculate T for each time
        # -------------------------
        # si z < zi, then we are in the solid (diffusion), else we are in the
        # liquid (T=cste).
        for t in range(1, len(self.time)):
            self.T[self.z < self.zi[t], t] = (
                (self.Tm - self.Tw) *
                erf(
                    self.z[self.z < self.zi[t]] / (
                        2*sqrt(self.alpha['solid']*self.time[t])
                    )
                ) /
                erf(self.delta) +
                self.Tw
            )
            self.T[self.z > self.zi[t], t] = (
                (self.Tm - self.Tinf) * erfc(
                    self.z[self.z > self.zi[t]] /
                    (2*sqrt(self.alpha['liquid']*self.time[t]))
                ) /
                erfc(
                    self.delta*sqrt(self.alpha['solid']/self.alpha['liquid'])
                ) +
                self.Tinf
            )
            self.T[self.z == self.zi[t], t] = self.Tm


class Ste_1d_1c_2p_qsqs_conf(object):
    def __init__(self):
        None

    def solve(self, Ta):
        time = linspace(0, 10, 100)

        mat = {'Water': ma.Water(), }
        mat = mat['Water']
        mat.set_temperature(mat.Tm)
        mat.set_phase('solid')

        cp = mat.get_cp()
        Lf = mat.get_Lf()
        ki = mat.get_k()
        Tm = mat.get_Tm()
        Ste = cp*(Tm - -15)/Lf
        rhoi = mat.get_rho()
        ai = mat.get_alpha()

        mat.set_phase('liquid')
        kl = mat.get_k()
        cpl = mat.get_cp()
        rhol = mat.get_rho()
        al = mat.get_alpha()

        aa = 20e-6
        ka = 0.024

        z0 = 1e-3

        zi_adim = [0]
        Tint = [2]
        Ta = 20

        time_adim = time*ai/z0**2

        dt = time_adim[1] - time_adim[0]

        # zi_adim.append(sqrt(2*Ste*dt))
        for n in range(1, len(time)):
            Stel = cpl*(Tint[n-1] - Tm)/Lf
            temp = zi_adim[n-1]**2 + 2*dt*(
                Ste - rhol/rhoi*al/ai*Stel/(1-zi_adim[n-1])*zi_adim[n-1]
            )

            zi_adim.append(sqrt(temp))
            # Tint_new = (
            #     Ta /
            #     (
            #         1 - kl/ka *
            #         erfc(z0/(2*sqrt(aa*time[n]))) * exp(z0**2/(4*aa*time[n])) *
            #         sqrt(pi*aa*time[n]) / (z0 - zi_adim[n-1]*z0)
            #     )
            # )

            Tint_new = (
                (
                    Tm +
                    Ta * ka/kl /
                    (sqrt(pi*aa*time[n])) *
                    exp(-z0**2/(4*aa*time[n])) /
                    erfc(z0/(2*sqrt(aa*time[n]))) *
                    (z0 - zi_adim[n]*z0)
                ) /
                (
                    1 - ka/kl * (z0 - zi_adim[n]*z0) /
                    sqrt(pi*aa*time[n]) *
                    exp(-z0**2/(4*aa*time[n])) /
                    erfc(z0/(2*sqrt(aa*time[n])))
                )
            )

            Tint.append(Tint_new)

        zi = array(zi_adim)*z0
        # zi = zi*(zi < z0) + z0*(zi >= z0)

        print('alpha rapport (l/i)', al/ai)
        print('rho rapport (l/i)', rhol/rhoi)
        print('Ste rapport (l/i)', Stel/Ste)
        print('rho*alpha rapport (l/i)', al/ai*rhol/rhoi)
        print('rho*alpha*Ste rapport (l/i)', al/ai*rhol/rhoi*Stel/Ste)

        # plt.plot(time, zi*1e3, '--k',)

        return zi, time, Tint


class Ste_1d_1c_2p_qsqs_conf0(object):
    def __init__(self):
        None

    def solve(self, Tint=20, Tw=-15):
        time = linspace(0, 10, 100)

        mat = {'Water': ma.Water(), }
        mat = mat['Water']
        mat.set_temperature(mat.Tm)
        mat.set_phase('solid')

        cp = mat.get_cp()
        Lf = mat.get_Lf()
        ki = mat.get_k()
        Tm = mat.get_Tm()
        Ste = cp*(Tm - Tw)/Lf
        rhoi = mat.get_rho()
        ai = mat.get_alpha()

        mat.set_phase('liquid')
        kl = mat.get_k()
        cpl = mat.get_cp()
        rhol = mat.get_rho()
        al = mat.get_alpha()
        Stel = cpl*(Tint - Tm)/Lf

        aa = 20e-6

        z0 = 1e-3

        zi_adim = [0]

        time_adim = time*ai/z0**2

        dt = time_adim[1] - time_adim[0]

        # zi_adim.append(sqrt(2*Ste*dt))
        for n in range(1, len(time)):
            temp = zi_adim[n-1]**2 + 2*dt*(
                Ste - rhol/rhoi*al/ai*Stel/(1-zi_adim[n-1])*zi_adim[n-1]
            )
            zi_adim.append(sqrt(temp))

        zi = array(zi_adim)*z0
        zi = zi*(zi < z0) + z0*(zi >= z0)

        print('alpha rapport (l/i)', al/ai)
        print('rho rapport (l/i)', rhol/rhoi)
        print('Ste rapport (l/i)', Stel/Ste)
        print('rho*alpha rapport (l/i)', al/ai*rhol/rhoi)
        print('rho*alpha*Ste rapport (l/i)', al/ai*rhol/rhoi*Stel/Ste)

        # plt.plot(time, zi*1e3, '--k',)

        return zi, time


class StefanMassique(object):
    """Stefan 1D monophasic model with fraction massique."""
    def __init__(self, mass_fraction):
        """Inite.

        add the mass fraction.
        """

        super().__init__()
        self.mass_fraction = mass_fraction

    def import_material(self, material):
        """
        Import Material. + add water cp to recalculate the initial temperature.
        """
        super().import_material(material)
        self.cp_l = ma.Water().get_cp()

    def import_geometry(self, z0, Nz, T_down, T_up, T_ini):
        super().import_geometry(z0, Nz, T_down, T_up, T_ini)
        self.T_down = T_up - self.mass_fraction*self.Lf/self.cp_l
        self.dTs = T_up - self.T_down

    def solve_delta(self):
        delta = sqrt(self.calc_St() / (2*(1-self.mass_fraction)))
        return delta


class PseudoTwoDim(object):
    def __init__(self):
        super().__init__()
        self.Lf = 333500

    def import_material(self, material):
        """Import the material properties."""
        if material is 'ice':
            mat = ma.Ice()
        elif material is 'water':
            mat = ma.Water()

        self.rho = mat.get_rho()
        self.k = mat.get_k()
        self.cp = mat.get_cp()
        self.alpha = mat.get_alpha()

    def import_geometry(self, z0, Nz, T_down, T_up, T_ini, r0):
        """Import geometry, boundaries condition and initiale condition.

        Parameters
        ----------
        z0 : float
            initiale height of the drop (in mm)
        Nz : int
            Number of node for z discretization
        T_down : float
            temperature in the substrate
        T_up : float
            temperature in the liquid
        T_ini : float
            initiale temperature
        rd : float
            Radius of the drop
        theta : float
            Contact angle of the drop (in degree)
        r0 : float
            Foot radius of the drop (in mm)

        """
        self.z0 = z0
        self.T_down = T_down
        self.dTs = T_up - T_down
        self.z = linspace(0, 1, Nz)
        self.T_ini = T_ini
        self.r0 = r0
        self.rd = (z0**2 + r0**2)/(2*z0)
        print('rd : ', self.rd)
        print('theta = ', arccos((self.rd - self.z0) / self.rd)*180/pi)
        self.theta = arccos((self.rd - self.z0) / self.rd)*180/pi

    def define_stop_condition(self, cdt='z0', dt=.1):
        """Redefine the stop condition."""
        self.dt = dt

    def front_position(self, adim=False):
        """Solve the front position."""

        """
        Euler solution
        """
        z_sq = [0]
        t = [0]
        a = 2*self.k*self.dTs / (self.rho*self.Lf)
        print('a = ', a)
        b = 1 / (self.r0**2)
        print('b = ', b)
        c = 2*b*(self.z0 - self.rd)
        print('c = ', c)
        it = 1
        singularity = c/(2*b)*(1 - sqrt(1+4*b/c**2))
        print('Singularity : ', singularity, 'z0 = ', self.z0)

        while sqrt(z_sq[it-1]) < self.z0 and t[it-1] < 300:
            t.append(it*self.dt)
            z_sq.append(
                z_sq[it-1] +
                a*self.dt /
                (1 - b*z_sq[it-1] + c*sqrt(z_sq[it-1]))
            )
            # print(sqrt(z_sq[it-1]), sqrt(z_sq[it]))
            it += 1

        z = sqrt(z_sq)
        self.time = t
        return z/self.z0 if adim else z

class PseudoTwoDim_trans(object):
    def __init__(self):
        super().__init__()
        self.Lf = 333500
        mat = ma.Ice()
        self.rho = mat.get_rho()
        self.k = mat.get_k()
        self.cp = mat.get_cp()
        self.alpha = mat.get_alpha()

        mat = ma.Water()
        self.rho_l = mat.get_rho()

    def set_geom(self, z0, r0, Nz):
        """Import geometry, boundaries condition and initiale condition.

        Parameters
        ----------
        z0 : float
            initiale height of the drop (in mm)
        Nz : int
            Number of node for z discretization
        r0 : float
            Foot radius of the drop (in mm)

        """
        self.z0 = z0
        self.z = linspace(0, z0, Nz)
        self.r0 = r0
        self.rd = (z0**2 + r0**2)/(2*z0)
        print('rd : ', self.rd)

    def set_therm(self, Tw, Tm):
        self.Tw = Tw
        self.Tm = Tm
        self.Ste = self.cp*(Tm - Tw)/self.Lf
        print(self.Ste)

    def set_time(self, t_end=30, Nt=1000):
        """Redefine the stop condition."""
        self.time = linspace(0, t_end, Nt)
        self.delta_t = self.time[1] - self.time[0]

    def get_zpos(self, adim=False):
        """Solve the front position."""

        """
        Euler solution
        """
        z = [0]
        delta = [sqrt(self.Ste/2)]
        n = 0

        b = 1 / (self.r0**2)
        c = 2*b*(self.z0 - self.rd)
        zf = self.rho_l/self.rho * self.z0

        while z[n] < zf:
            z1 = sqrt(
                z[n]**2 +
                4*self.delta_t * self.k * (self.Tm - self.Tw) / (self.rho*self.Lf) *  # diffusion therm
                exp(-delta[n]**2) / erf(delta[n]) * delta[n] /  # transitoire correction
                (1 - b*z[n]**2 + c*z[n])  # geometric correction
            )

            delta1 = z1/(2*sqrt(self.alpha*(n+1)*self.delta_t))

            z.append(z1)
            delta.append(delta1)
            n += 1

        # def get_delta(zn):
        #     b = 1 / (self.r0**2)
        #     c = 2*b*(self.z0 - self.rd)
        #     right_delta = self.Ste/sqrt(pi)*(
        #         1-b*zn*2+c*zn
        #     )**(-1)
        #
        #     def func(delta):
        #         return delta*exp(delta**2)*erf(delta) - right_delta
        #
        #     return op.fsolve(func, .1)
        #
        # while n < len(self.time)-1:
        #     n += 1
        #     print('delta : ', get_delta(z[n-1])[0])
        #     delta.append(get_delta(z[n-1]))
        #     z.append(2*delta[0]*sqrt(self.alpha*self.time[n]))

        return z/self.z0 if adim else z


if __name__ == '__main__':
    # case = {
    #     0: Ste_1d_1c_1p_qs, 1: Ste_1d_1c_1p_is,
    #     2: Ste_1d_1c_2p_qsis, 3: Ste_1d_1c_2p_isis,
    #     4: Ste_1d_1c_2p_isad
    # }
    # colormaps = {
    #     0: 'Blues', 1: 'Reds', 2: 'Oranges', 3: 'Greens', 4: plt.cm.bone
    # }
    # colors = {
    #     0: 'tab:blue', 1: 'tab:red', 2: 'tab:orange', 3: 'tab:green', 4: 'k'
    # }
    #
    # fig1 = plt.figure(figsize=[8, 4.5])
    # ax0 = plt.subplot(1, 1, 1)
    # ax0.grid(True)
    # ax0.set_xlabel('t (s)')
    # ax0.set_ylabel(r'$z_i$ (mm)')
    #
    # fig3 = plt.figure(figsize=[8, 4.5])
    # ax30 = plt.subplot(1, 3, 1)
    # ax30.grid(True)
    # ax30.set_xlabel('T (°C)')
    # ax30.set_ylabel(r'$z$ (mm)')
    # ax30.set_title('Initiale time')
    #
    # ax31 = plt.subplot(1, 3, 2)
    # ax31.grid(True)
    # ax31.set_yticklabels([])
    # ax31.set_xlabel('T (°C)')
    # ax31.set_title('Half time')
    #
    # ax32 = plt.subplot(1, 3, 3)
    # ax32.set_yticklabels([])
    # ax32.grid(True)
    # ax32.set_xlabel('T (°C)')
    # ax32.set_title('End time')
    #
    # dil = False
    # for i in [0, 0, 1, 1]:
    #     ste = case[i](dilatation=dil)
    #     dil = not dil
    #
    #     ste.set_material()
    #     ste.set_geometry(z0=2e-3, zmax=3e-3, Nz=200)
    #     ste.set_thermic(T=-15)
    #     ste.set_time(tmax=30, Nt=1000)
    #     ste.solve()
    #     T = ste.get_tempfield()
    #     zi = ste.get_front()
    #     zt = ste.get_total_height()
    #     zf = ste.get_zf()
    #     tf = ste.get_tf()
    #     z, t = ste.get_mesh(), ste.get_time()
    #     print(
    #         i, ste.get_Ste(),
    #         ste.get_tz0(), ste.get_tf(),
    #         ste.get_delta(),
    #         -15, None
    #     )
    #     ax0.plot(t, zi*1e3, '-', c=colors[i])
    #     ax0.plot([tf, tf], [0, zf*1e3], '.--', c=colors[i])
    #     ax0.plot([0, tf], [zf*1e3, zf*1e3], '.--', c=colors[i])
    #
    #     fig2 = plt.figure(figsize=[8, 4.5])
    #     ax1 = plt.subplot(1, 1, 1)
    #     cs = ax1.contourf(
    #         t, z*1e3, T,
    #         cmap=colormaps[i], levels=[-15, -10, -5, -1e-3, 0]
    #     )
    #     ax1.plot(t, zt*1e3, '--', c='k')
    #     ax1.plot(t, zi*1e3, '--', c=colors[i])
    #     fig2.colorbar(cs)
    #
    #     ax30.plot(T[:, 1], z, ls='-.', c=colors[i])
    #     ax31.plot(T[:, int(len(T[:, 0])/2)], z, ls='-.', c=colors[i])
    #     ax32.plot(T[:, -1], z, ls='-.', c=colors[i])
    #
    # for i in [2, 2, 3, 3, 4, 4]:
    #     ste = case[i](dil)
    #     dil = not dil
    #
    #     ste.set_material()
    #     ste.set_geometry(z0=2e-3, zmax=3e-3, Nz=200)
    #     ste.set_thermic(T=[-15, 10])
    #     ste.set_time(tmax=30, Nt=1000)
    #     ste.solve()
    #     T = ste.get_tempfield()
    #     zi = ste.get_front()
    #     zt = ste.get_total_height()
    #     z, t = ste.get_mesh(), ste.get_time()
    #
    #     print(
    #         i, ste.get_Ste(),
    #         ste.get_tz0(), ste.get_tf(),
    #         ste.get_delta(),
    #         -15, 10
    #     )
    #
    #     ax0.plot(t, zi*1e3, '-', c=colors[i])
    #     ax0.plot([ste.tz0, ste.tz0], [0, ste.z0*1e3], '.--', c=colors[i])
    #     ax0.plot([0, ste.tz0], [ste.z0*1e3, ste.z0*1e3], '.--', c=colors[i])
    #
    #     fig2 = plt.figure(figsize=[8, 4.5])
    #     ax1 = plt.subplot(1, 1, 1)
    #     cs = ax1.contourf(
    #         t, z*1e3, T,
    #         cmap=colormaps[i], levels=[-15, -10, -5, 0, 5, 10]
    #     )
    #     ax1.plot(t, zt*1e3, '--', c='k')
    #     ax1.plot(t, zi*1e3, '--', c='tab:orange')
    #     fig2.colorbar(cs)
    #
    #     ax30.plot(T[:, 1], z, ls='-.', c=colors[i])
    #     ax31.plot(T[:, int(len(T[:, 0])/2)], z, ls='-.', c=colors[i])
    #     ax32.plot(T[:, -1], z, ls='-.', c=colors[i])
    #
    # plt.show()
    #
    # ste = case[3](dilatation=False)
    # ste.set_material()
    # ste.set_geometry(z0=2e-3, zmax=3e-3, Nz=200)
    # ste.set_thermic(T=[-15, 10])
    # ste.set_time(tmax=30, Nt=1000)
    # ste.solve()
    # ste.plot_front()
    # ste.plot_heat_field()

    ste = Ste_1d_1c_1p_qs(dilatation=False)
    ste.set_material()
    ste.set_geometry(z0=1e-3, zmax=3e-3, Nz=200)
    ste.set_thermic(T=-15)
    ste.set_time(tmax=10, Nt=1000)
    ste.solve()
    T = ste.get_tempfield()
    zi = ste.get_front()
    zt = ste.get_total_height()
    zf = ste.get_zf()
    tf = ste.get_tf()
    z, t = ste.get_mesh(), ste.get_time()
    plt.plot(t, zi*1e3, '-', c='red', label='Analytical solution')
    # tau = 3.5
    # t = linspace(0, 20, 100)
    # plt.plot(t, 1/4*sqrt(t)+1/2*(1-exp(-t/tau)**2), '--', color='tab:green')
    ste = Ste_1d_1c_2p_qsqs_conf0()
    colors = ['k', 'tab:purple', 'tab:blue', 'tab:green', 'tab:orange', 'tab:red']
    phis = ['0%', '30%', '40%', '50%', '60%', '70%']
    i = 0
    for Tint in [0, 2, 6, 9, 12, 14]:
        zi, time = ste.solve(Tint, -15)
        plt.plot(
            time, zi*1e3, '-.',
            color=colors[i], label=r'$\phi$ : ' + phis[i]
        )
        i += 1

    ste = Ste_1d_1c_2p_qsqs_conf()
    zi, time, Tint = ste.solve(20)

    plt.plot(time, zi*1e3, '-', color='tab:red')

    plt.grid(True)
    plt.legend(fancybox=True, shadow=True, loc='lower right')
    plt.xlabel('t (s)')
    plt.xlabel('z (mm)')

    plt.figure()
    plt.plot(time, Tint, '-.k')
    plt.grid(True)

    plt.show()
