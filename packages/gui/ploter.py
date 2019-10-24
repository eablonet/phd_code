# redefine plot to create beautiful figure()

from matplotlib import pyplot as plt
from matplotlib import ticker as tck


class Figure:
    """Redefine figure."""

    def __init__(self):
        """Create a figure."""
        self.fig = plt.figure(figsize=(8, 4.5), constrained_layout=True)
        self.ax = []

    def __len__(self):
        """Return the number of axes."""
        return len(self.ax)

    def set_win_title(self, txt):
        """Set the title of the figure."""
        self.fig.canvas.set_window_title(txt)

    def show(self):
        """Display the figure."""
        plt.show()

    def add_ax(self, n):
        """Convert add_subplot method."""
        self.ax.append(Ax(self.fig.add_subplot(n)))
        return self.ax[-1]

    def colorbar(self, mapname, values, clabel=None, ticks=None, ticklabels=None):
        """Add a colorbar."""
        cmap = plt.cm.get_cmap(mapname, len(values))
        # Make dummie mappable
        x, _ = self.ax[-1].ax.get_xlim()
        y, _ = self.ax[-1].ax.get_ylim()
        x = [x for _ in range(len(values))]
        y = [y for _ in range(len(values))]
        dummie_cax = self.ax[-1].ax.scatter(
            x, y, c=range(len(values)), cmap=cmap
        )

        cbar = self.fig.colorbar(dummie_cax)
        dummie_cax.set_visible(False)

        if clabel is not None:
            cbar.set_label(clabel)

        cbar.set_ticks([0, len(values)-1])
        cbar.set_ticklabels([values[0], values[-1]])


class Ax:
    """Redefine ax properties."""

    def __init__(self, ax):
        """Create the ax."""
        self.ax = ax
        self.ax.grid(True)

        self._spines()

    def _spines(self):
        """Redesign spines."""
        self.ax.spines['right'].set_visible(False)
        self.ax.spines['top'].set_visible(False)

    def plot(self, *args, **kwargs):
        """Add a plot."""
        self.ax.plot(*args, lw=2, ms=8, **kwargs)

    def loglog(self, *args, **kwargs):
        """Add a loglog."""
        self.ax.loglog(*args, lw=2, ms=8, **kwargs)

    def legend(self, **kwargs):
        """Enable legend."""
        self.ax.legend(
            fancybox=True, shadow=True,
            **kwargs
        )

    def imshow(self, *args, **kwargs):
        """Add a loglog."""
        self.ax.imshow(*args, cmap='gray', **kwargs)

    def xlabel(self, txt):
        """Set xlabel."""
        self.ax.set_xlabel(txt, fontsize=14)

    def ylabel(self, txt):
        """Set ylabel."""
        self.ax.set_ylabel(txt, fontsize=14)

    def title(self, *args, **kwargs):
        """Define axe title."""
        self.ax.set_title(*args, fontsize=14, **kwargs)


if __name__ == '__main__':
    fig = Figure()
    fig.set_win_title('Gallery figure')

    # data
    x = [i/100 for i in range(-50, 150)]
    y = [i**3 for i in x]

    # create ax
    ax = fig.add_ax(111)

    # plot
    ax.plot(x, y, ls='--', marker='o', color='tab:blue', markevery=10, label=r'$y^3$')
    ax.plot(x, x, ls='-', color='tab:red', marker='v', markevery=20, label=r'$y$')

    ax.ax.set_title("Comparison of power")

    # add legend
    ax.legend(loc='upper left')

    # label
    ax.xlabel('Position')
    ax.ylabel('Energie')

    # display figure
    fig.show()

    # colorbar
    fig = Figure()
    fig.set_win_title('Gallery figure')
    ax = fig.add_ax(111)
    import numpy as np
    t = np.linspace(0, 2, 100)
    cmap = plt.get_cmap('tab10', 4)
    for T in range(1, 5):
        y = np.cos(2*np.pi*t/T)
        ax.plot(
            t, y,
            color=cmap(T-1)
        )
    fig.colorbar(
        mapname='tab10', values=range(1, 5),
        ticks=range(1, 5),
        ticklabels=['1']
    )
    fig.show()
