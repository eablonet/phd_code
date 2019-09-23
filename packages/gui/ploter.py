# redefine plot to create beautiful figure()

from matplotlib import pyplot as plt
from matplotlib import ticker as tck


class Figure:
    """Redefine figure."""

    def __init__(self):
        """Create a figure."""
        self.fig = plt.figure(figsize=(8, 4.5))
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

    def legend(self, **kwargs):
        """Enable legend."""
        self.ax.legend(
            fancybox=True, shadow=True,
            **kwargs
        )

    def xlabel(self, txt):
        """Set xlabel."""
        self.ax.set_xlabel(txt, fontsize=14)

    def ylabel(self, txt):
        """Set ylabel."""
        self.ax.set_ylabel(txt, fontsize=14)


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
