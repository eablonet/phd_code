"""
Stefan Solver Class.

----
The main class is the generic Stefan class.
Other class are depending on condition.
----
class
  Stefan : main

"""

from numpy import sqrt, linspace, zeros, pi, exp, arccos, arange
from numpy import cos
from matplotlib import pyplot as plt
from scipy.special import erf, erfc
from scipy import optimize as op

try:
    from packages.main import Material as ma
except ModuleNotFoundError:
    import Material as ma


class Stefan(object):
    """Main class.

    This class provide the One Dimensional steady monophasic solution of Stefan
    problem.

    Parameters
    ----------
    rho :
        material density
    k :
        thermal conductivity
    cp :
        specific capacity
    alpha :
        diffusivity

    Methods
    -------
    import_material:
        import value of rho, k, cp and alpha

    """

    def __init__(self):
        """Allocate variables."""
        self.rho = None
        self.k = None
        self.cp = None
        self.alpha = None
        self.Lf = 333500

        self.dTs = None
        self.T_down = None
        self.T_ini = None

        self.time = None
        self.z = None

        self.z0 = None

    def __repr__(self):
        """Print function."""
        return (
            'One dimensional steady monophasic solution of Stefan Problem.'
        )

    def import_material(self, material):
        """Import the material properties."""
        if material is 'ice':
            mat = ma.Water('solid')
        elif material is 'water':
            mat = ma.Water('liquid')

        self.rho = mat.get_rho()
        self.k = mat.get_k()
        self.cp = mat.get_cp()
        self.alpha = mat.get_alpha()

    def import_geometry(self, z0, Nz, T_down, T_up, T_ini):
        """Import geometry, boundaries condition and initiale condition."""
        self.z0 = z0
        self.T_down = T_down
        self.dTs = T_up - T_down
        self.z = linspace(0, 1, Nz)
        self.T_ini = T_ini

    def calc_St(self):
        """Determine Stefan number with variable."""
        return self.cp*self.dTs/self.Lf

    def solve_delta(self):
        """Calculate delta value.

        Whatever is adimentional or not delta have the same value.
        """
        return sqrt(self.calc_St()/2)

    def calc_tf(self, adim=False):
        """Calculate the tf.

        Parameters
        ----------
        adim : bool
            if True, return the adimentional temperature,
            if False, return the dimensional one.

        """
        tf = 1 / (4*self.solve_delta()**2)
        return tf if adim else tf * self.z0**2/self.alpha

    def calc_t12(self, adim=False):
        """Calculate the half-life time."""
        if adim:
            return 1/(16*self.solve_delta()**2)
        else:
            t0 = self.z0**2/self.alpha
            return t0/(16*self.solve_delta()**2)

    def define_stop_condition(self, cdt='tf', Nt=100, val=1e-3):
        """Choose the end condition.

        In function of the method chosen, determine the stop time for
        calculation.

        Parameters
        ----------
        cdt : string
            Condition d'arrêt pour le givrage.
            A choisir parmi : tf, z0, zf, zp:
                tf : theoritical final time
                z0 : when z is egual to the initiale height
                zf : when z is egual to theoritical final height (equivalent to
                    tf)
                zp : personal height value (in m)
                tp : personal time
        Nt : int
            Number of steps for discretize the time.
        val : float
            Value for zp or tp if one of this conditions is used.

        """
        if cdt in ['tf', 'zf']:
            tf = self.calc_tf()
        elif cdt == 'z0':
            tf = 1/(4*self.solve_delta()**2)
        elif cdt == 'zp':
            tf = val**2/self.z0**2/(4*self.solve_delta()**2)  # val == zp
        elif cdt == 'tp':
            tf = val  # val == tp

        self.time = linspace(0, tf, Nt)

    def front_position(self, adim=False):
        """Calculate front position over time.

        Parameters
        ----------
        adim : bool
            if True, return the adimentional temperature,
            if False, return the dimensional one.

        """
        z = 2*self.solve_delta()*sqrt(self.alpha*self.time)
        if adim:
            z /= self.z0

        return z

    def temp_field(self, adim=False):
        """Calculate the temperature field.

        Give the temperature T(z, t) at each position of the domain
        for each time.

        Parameters
        ----------
        adim : bool
            if True, return the adimentional temperature,
            if False, return the dimensional one.

        """
        T = zeros([len(self.z), len(self.time)])
        T[:, 0] = (self.T_ini - self.T_down) / self.dTs
        zs_star = self.front_position(adim=True)

        for t in range(1, len(self.time)):
            for z in range(len(self.z)):
                if self.z[z] < zs_star[t]:
                    T[z, t] = self.z[z]/zs_star[t]
                else:
                    T[z, t] = 1

        return T if adim else T*self.dTs + self.T_down

class StefanMassique(Stefan):
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

class Diphasique(Stefan):
    """Solution diphasique 1D analytique avec dilatation."""

    def __init__(self):
        """Re-Allocate the variables."""
        super().__init__()

        # new variables
        self.dTl = None
        self.T_m = None

        self.rho_l = None
        self.cp_l = None
        self.k_l = None
        self.alpha_l = None

    def __repr__(self):
        """Print function."""
        return (
            'One dimensional diphasic unsteady solution of Stefan problem ' +
            'with dilatation'
        )

    def import_geometry(self, z0, Nz, T_down, T_up, T_ini, T_m):
        """Redefine the geometry function."""
        super().import_geometry(z0, Nz, T_down, T_m, T_ini)
        self.dTl = T_m - T_up
        self.T_m = T_m

    def import_material(self, mat1, mat2):
        """Surcharge de import material pour 2 materiaux."""
        super().import_material(mat1)

        if mat2 is 'ice':
            mat = ma.Ice()
        elif mat2 is 'water':
            mat = ma.Water()

        self.rho_l = mat.get_rho()
        self.k_l = mat.get_k()
        self.cp_l = mat.get_cp()
        self.alpha_l = mat.get_alpha()

    def calc_St_liquid(self):
        """Calculate the liquid Stefan number."""
        return self.cp_l*self.dTl / self.Lf

    def calc_tf(self, adim=False):
        """Calculate the tf.

        Parameters
        ----------
        adim : bool
            if True, return the adimentional temperature,
            if False, return the dimensional one.

        """
        tf = (self.rho_l/self.rho)**2/(4*self.solve_delta()**2)
        return tf if adim else tf*self.z0**2/self.alpha

    def solve_delta(self):
        """Solve the value of delta."""
        a = self.rho/self.rho_l * sqrt(self.alpha/self.alpha_l)

        def func(delta):
            f = (
                sqrt(pi)*delta*exp(delta**2)*erf(delta) -
                self.calc_St() -
                self.calc_St_liquid()*erf(delta)*exp(delta**2*(1-a**2)) / (
                    a*(erfc(a*delta))
                )
            )
            return f

        return op.fsolve(func, .1)

    def temp_field(self, adim=False):
        """Calculate the temperature field.

        Give the temperature T(z, t) at each position of the domain
        for each time.

        Parameters
        ----------
        adim : bool
            if True, return the adimentional temperature,
            if False, return the dimensional one.

        """
        Ts = zeros([len(self.z), len(self.time)])
        Tl = zeros([len(self.z), len(self.time)])
        Tl[:, 0] = (self.T_ini - self.T_m) / self.dTl
        # Ts[:, 0] = (self.T_ini - self.T_down) / self.dTs
        zs_star = self.front_position(adim=True)
        delta = self.solve_delta()

        for t in range(1, len(self.time)):
            for z in range(len(self.z)):
                eta = self.z[z] / (2*sqrt(self.time[t]))
                if self.z[z] < zs_star[t]:
                    # in the first material
                    Ts[z, t] = erf(eta)/erf(delta)
                    # Tl[z, t] = 0
                elif self.z[z] > zs_star[t]:
                    # in the second material
                    # Ts[z, t] = 0
                    Tl[z, t] = (
                        erfc(
                            sqrt(self.alpha/self.alpha_l)*(
                                eta - (1-self.rho/self.rho_l)*delta
                            )
                        ) /
                        erfc(
                            sqrt(self.alpha/self.alpha_l) * delta *
                            self.rho/self.rho_l
                        )
                    )
                    # print('delta = ', delta, ' eta = ', eta, ' Temp = ', Tl[z, t])
                else:
                    # in the front line
                    # Tl[z, t] = 0
                    Ts[z, t] = 1

        return Ts+Tl if adim else Ts*self.dTs + self.T_down + self.dTl + \
            self.T_m


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
    st = Stefan()
    st.import_material('ice')
    st.import_geometry(z0=2e-3, Nz=100, T_down=-15, T_up=0, T_ini=10)
    st.define_stop_condition('tf', Nt=1000)

    # st_m = StefanMassique(.15)
    # st_m.import_material('ice')
    # st_m.import_geometry(z0=2e-3, Nz=100, T_down=-10, T_up=0, T_ini=10)
    # print('Nucleation temp correction', st_m.T_down)
    # st_m.define_stop_condition('tf', Nt=1000)
    #
    # st_m2 = StefanMassique(.05)
    # st_m2.import_material('ice')
    # st_m2.import_geometry(z0=2e-3, Nz=100, T_down=-10, T_up=0, T_ini=10)
    # print('Nucleation temp correction', st_m2.T_down)
    # st_m2.define_stop_condition('tf', Nt=1000)

    # dip = Diphasique()
    # dip.import_geometry(z0=2e-3, Nz=100, T_down=-15, T_up=10, T_ini=10, T_m=0)
    # dip.import_material('ice', 'water')
    # dip.define_stop_condition('tf')

    ptd = PseudoTwoDim()
    ptd.import_geometry(z0=2e-3, Nz=100, T_down=-15, T_up=0, T_ini=10, r0=2e-3/.64)
    ptd.import_material('ice')
    ptd.define_stop_condition('tf', dt=.01)
    z = ptd.front_position()  # to generate the time

    # ptd_new = PseudoTwoDim_trans()
    # ptd_new.set_geom(2e-3, 2e-3/.64, 1000)
    # ptd_new.set_therm(Tw=-15, Tm=0)
    # ptd_new.set_time(t_end=30, Nt=1000)
    # ptd_new.set_time(t_end=30, Nt=1000)
    # z_new = ptd_new.get_zpos()

    plt.figure(figsize=[8, 4.5])
    zs_comp = sqrt(2*2.2*(0+15)/(916.2*333500)*st.time)
    # plt.plot(st.time, zs_comp, '--m', label='Ref')
    plt.plot(st.time, st.front_position(), '--k', label='Quasi-statique')
    # plt.plot(dip.time, dip.front_position(), '--r', label='Full Stefan')
    plt.plot(ptd.time, z, '--', c='tab:orange', label='Geom a = .64')
    # plt.plot(arange(len(z_new))*ptd_new.delta_t, z_new, '--', c='tab:orange', label='Pseudo-2D_trans')
    # plt.plot(st_m.time, st_m.front_position(), '--g', label='Stefan massique fraction .15')
    # plt.plot(st_m2.time, st_m2.front_position(), '--c', label='Stefan massique fraction .3')

    ptd = PseudoTwoDim()
    ptd.import_geometry(z0=2e-3, Nz=100, T_down=-15, T_up=0, T_ini=10, r0=2.01e-3)
    ptd.import_material('ice')
    ptd.define_stop_condition('tf', dt=.01)
    z = ptd.front_position()  # to generate the time
    plt.plot(ptd.time, z, '--', c='tab:blue', label='Geom hemi-sphérique')

    plt.legend(fancybox=True)
    plt.grid(True)
    plt.xlabel('t (s)')
    plt.ylabel(r'$z_s$ (mm)')

    plt.figure()
    plt.plot(st.z, st.temp_field(True)[:, 30])
    # plt.plot(dip.z, dip.temp_field(True)[:, 30])
    # plt.plot(
    #     [dip.front_position(adim=True)[30], dip.front_position(adim=True)[30]],
    #     [0, 1],
    #     '--b'
    # )

    plt.figure()
    plt.plot(st.time, st.temp_field(True)[30, :])
    # plt.plot(dip.time, dip.temp_field(True)[30, :])
    plt.show()
