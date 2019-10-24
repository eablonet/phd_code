# diffusion solver

class Mesh:
    def __init__(self):
        # z-space
        self.dz = None
        self.H = None
        self.z = None

        # t-space
        self.dt = None
        self.t = None
        self.ts = None
        self.tend = None

    @property
    def Nz(self):
        """Return the number of point for z discretization."""
        return int(self.H[1]/self.dz)

    @property
    def Nt(self):
        """Return the number of point for t discretization."""
        return int((self.tend - self.ts)/self.dt)

    def get_mesh(self):
        """Generate mesh."""
        self.Nz = int(self.H[1]/self.dz)
        self.z = arange(0, self.H[1], self.dz[1])
        return self.z

    def get_time_discretization(self):
        """Generate time discretization."""
        self.t = arange(self.ts, self.tend, self.dt)
        return self.t

class Field:
    """Class for thermal field."""
    def __init__(self):
        """Allocate vectors."""
        self.mesh = None
        self.zi = zeros(self.mesh.Nt)  # front location
        self.alpha = zeros((self.mesh.Nz, self.mesh.Nt))  # phase function (0: liquid; 1: ice)
        self.T = zeros((self.mesh.Nz, self.mesh.Nt))

    def add_thermal(self, T, j, n):
        """Add point."""
        self.T[j, n] = T

    def add_front_point(self, zi, n):
        self.zi[n] = zi

    def add_alpha_point(self, alpha, j, n):
        self.alpha[j, n] = alpha

    def plot(self):
        """Display ThermalField."""
        plt.contourf(self.t, self.z, self.T)
        plt.show()

class Diffusion:
    def __init__(self):
        self.mesh = Mesh()
        self.field = Field()
        None

    def time_step(self):
        self.new_front()

    def new_front(self, k, n):
        znew = self.field.zi[n-1] + \
            self.mesh.dt/(self.rho[2]*self.Lf*self.mesh.dz) * (
                self.field.T[k, n-1] * (self.k[1]+self.k[2]) -
                
            )
        return
