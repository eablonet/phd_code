# -*- coding: utf-8 -*-
"""
Create spline to track object.

Authors
-------
eablonet

"""

import numpy as np
from scipy import interpolate

from scipy.interpolate import UnivariateSpline

from matplotlib import pyplot as plt


class SplineBuilder(object):
    """
    Create Line on graph.

    Methods
    -------
    __call__
        on call the method plot the line

    """

    def __init__(
        self,
        ax1, ax2, ax3, ax4, image, grad, grady,
        ax, intensity, z_intensity, r_intensity
    ):
        """Do the initiation."""
        self.axes = [ax1, ax2, ax3, ax4]
        self.axes2 = ax

        self.image = image
        self.grad = grad
        self.grady = grady

        self.sel_ax_zoom = 0

        self.intensity = intensity
        self.z_intensity = z_intensity
        self.r_intensity = r_intensity

        self.press = False  # init press indicator

        self.dep_mode = False  # at init no point is selected
        self.selected = [False, []]
        # first indicate if there is a selected point,
        # if True the second indicata is index

        # line 1 #
        self.ax1 = ax1
        self.line1 = ax1.plot([], [])[0]
        self.line1.set_linestyle('-')
        self.line1.set_color('b')
        self.line1.set_linewidth(.8)

        self.pt1 = ax1.plot([], [])[0]
        self.pt1.set_linestyle('None')
        self.pt1.set_marker('o')
        self.pt1.set_markeredgecolor('b')

        self.dp1 = ax1.plot([], [])[0]
        self.dp1.set_linestyle('None')
        self.dp1.set_marker('s')
        self.dp1.set_markeredgecolor('r')
        self.dp1.set_markerfacecolor('none')

        # line 2 #
        self.line2 = ax2.plot([], [])[0]
        self.line2.set_linestyle('-')
        self.line2.set_color('b')
        self.line2.set_linewidth(.8)

        self.pt2 = ax2.plot([], [])[0]
        self.pt2.set_linestyle('None')
        self.pt2.set_marker('o')
        self.pt2.set_markeredgecolor('b')

        # line 3 #
        self.cross_x = ax3.plot(
            [], [], ls='-', color='tab:red', lw=2
        )[0]  # x(cste)
        self.cross_y = ax3.plot(
            [], [], ls='-', color='tab:red', lw=2
        )[0]  # y(cste)

        # line 4 #
        self.line4 = ax4.plot([], [])[0]
        self.line4.set_linestyle('-')
        self.line4.set_color('b')
        self.line4.set_linewidth(.8)

        self.pt4 = ax4.plot([], [])[0]
        self.pt4.set_linestyle('None')
        self.pt4.set_marker('o')
        self.pt4.set_markeredgecolor('b')

        self.dots = []
        for a in self.axes2:
            self.dots.append(
                a.plot(
                    [], [],
                    ls='None',
                    marker='o', ms=6, color='tab:red', mfc='None'
                )[0]
            )

        # enable connections
        # ------------------
        self.cid1 = self.ax1.figure.canvas.mpl_connect(
            'button_press_event', self.enable_press
        )
        self.cid1 = self.ax1.figure.canvas.mpl_connect(
            'button_release_event', self.on_press
        )
        self.cid1_key = self.ax1.figure.canvas.mpl_connect(
            'key_release_event', self.on_press_key
        )

        self.cid1_move = self.ax1.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion
        )

        self.cid2 = self.axes2[0].figure.canvas.mpl_connect(
            'button_release_event', self.on_press_intensity
        )

        self.xs = np.array(self.line1.get_xdata(), dtype=np.int)
        self.ys = np.array(self.line1.get_ydata(), dtype=np.int)
        self.xs_interp = np.array(self.line1.get_xdata(), dtype=np.int)
        self.ys_interp = np.array(self.line1.get_ydata(), dtype=np.int)

    def enable_press(self, event):
        self.press = True

    def update_fig(self):
        """Update the figure."""
        self.pt1.set_data(self.xs, self.ys)
        self.pt2.set_data(self.xs, self.ys)
        self.pt4.set_data(self.xs, self.ys)

        self.line1.set_data(self.xs_interp, self.ys_interp)
        self.line2.set_data(self.xs_interp, self.ys_interp)
        self.line4.set_data(self.xs_interp, self.ys_interp)

        self.line4.figure.canvas.draw()

    def interpolate(self):
        """Construct the line interpolation."""
        if len(self.xs) > 1:
            f = interpolate.interp1d(
                self.xs_interp, self.ys_interp,
                'slinear'
            )
            self.xs_interp = np.arange(
                np.min(self.xs), np.max(self.xs), 1
            )
            self.ys_interp = np.array(
                f(self.xs_interp), dtype=np.int
            )

    def on_press(self, event):
        """Create line on clicking on the graph.

        Actions
        -------
            Press left button to add a point,
            Middle button to finish the plot,
            Right to remove the last one.

        """
        print(event.button)
        if self.press is False:
            return
        else:
            self.press = False

        zoom_cond = (
            self.line1.axes.get_navigate_mode() == 'ZOOM' or
            self.line2.axes.get_navigate_mode() == 'ZOOM' or
            self.line4.axes.get_navigate_mode() == 'ZOOM'
        )

        if zoom_cond:
            return

        if event.inaxes not in [
            self.line1.axes, self.line2.axes, self.line4.axes,
        ]:
            return

        elif event.button == 2:
            """
            Clic du milieu pour valider les points et
            passer à l'image suivante.
            """
            plt.close()

        elif (
            event.button == 3 and
            len(self.xs) > 0 and len(self.ys) > 0
        ):
            """
            Clic droit supprime le dernier point.
            """
            if len(self.xs) > 0 and len(self.ys) > 0:
                self.xs = np.delete(self.xs, -1)
                self.ys = np.delete(self.ys, -1)

        elif event.button == 1:
            """
            Clic gauche, ajoute un point.
            Avec ce mode d'ajout, le code va retrier les points dans
            l'ordre des x croissant.
            Il ne sont donc pas prient en compte dans l'ordre d'apparition.
            """
            self.xs = np.append(self.xs, int(event.xdata))
            self.ys = np.append(self.ys, int(event.ydata))

        self.ys_interp = self.ys[np.argsort(self.xs)]
        self.xs_interp = np.sort(self.xs)

        self.interpolate()
        self.update_fig()

        # new quicker method
        if len(self.xs_interp) > 1:
            x_new = self.xs[-1]  # get last x created
            x_sort = np.sort(self.xs)  # sort xs
            if x_new == min(self.xs_interp):
                x_min = x_sort[0]
                x_max = x_sort[1]
            elif x_new == max(self.xs_interp):
                x_min = x_sort[-2]
                x_max = x_sort[-1]
            else:
                x_min = x_sort[(x_sort == x_new)-1]
                x_max = x_sort[(x_sort == x_new)+1]

            x_sort = self.xs_interp[x_min:x_max]
            for n in range(7):
                if self.r_intensity[n] in x_sort:
                    z = self.z_intensity[x_sort == self.r_intensity]
                    intens = self.intensity[z, n]
                    self.axes2[n].plot(
                        intens,
                        z, 'or', ms=6
                    )
                    self.axes2[n].figure.canvas.draw()

    def on_press_intensity(self, event):
        for a in self.axes2:
            if a.get_navigate_mode() == 'ZOOM':  # not zoom mode
                return

        if event.inaxes not in self.axes2:  # not in the ax
            return

        if event.button == 1:  # left clic
            x_curs = int(event.xdata)
            y_curs = int(event.ydata)

            # get axe number
            n = 0
            while event.inaxes != self.axes2[n]:
                n += 1

            # look for the closet point from clic
            d = np.inf  # init distance pt from line
            z_out = 0
            for z in self.z_intensity[y_curs-15:y_curs+15]:
                dx = x_curs - self.intensity[z, n]
                dy = y_curs - z
                if np.sqrt(dx**2 + dy**2) < d:
                    d = np.sqrt(dx**2 + dy**2)
                    z_out = z

            # find the local maximum [+- 4px]
            it_max = 0  # init intensity maximum
            z_max = 0  # init z of maximum intensity
            for z in self.z_intensity[z_out-4:z_out+4]:
                if self.intensity[z, n] > it_max:
                    it_max = self.intensity[z, n]
                    z_max = z

            self.dots[n].set_xdata(self.intensity[z_max, n])
            self.dots[n].set_ydata(z_max)
            self.axes2[n].figure.canvas.draw()

            # update fig
            self.xs = np.append(self.xs, int(self.r_intensity[n]))
            self.ys = np.append(self.ys, int(z_out))
            self.ys_interp = self.ys[np.argsort(self.xs)]
            self.xs_interp = np.sort(self.xs)

            self.interpolate()
            self.update_fig()

    def on_motion(self, event):
        if self.press is False:
            return
        if event.inaxes not in self.axes:
            return

        x1, x2 = int(event.xdata), int(event.xdata)
        y1, y2 = int(event.ydata), int(event.ydata)
        le = 20  # len of the zoom

        x1 = (x1-le)*(x1 > le)   # = x1-le si x1 > le 0 sinon 0
        y1 = (y1-le)*(y1 > le)

        x2 = (x2+le)*(x2 < np.shape(self.image)[1] - le+1) + np.shape(self.image)[1] * (1 - (x2 < np.shape(self.image)[1] - le+1))
        y2 = (y2+le)*(y2 < np.shape(self.image)[0] - le+1) + np.shape(self.image)[0] * (1 - (y2 < np.shape(self.image)[0] - le+1))

        zoom = (
            self.image * (self.sel_ax_zoom == 0) +
            self.grad * (self.sel_ax_zoom == 1) +
            self.grady * (self.sel_ax_zoom == 3)
        )  # select which image to use (change pressing 'o, 'g or 'y')

        zoom = zoom[y1:y2, x1:x2]  # crop the image

        xc = int(np.shape(zoom)[1]/2)  #  x_center_zoom
        yc = int(np.shape(zoom)[0]/2)  #  y_center_zoom

        self.cross_x.set_xdata([xc, xc])
        self.cross_y.set_ydata([yc, yc])

        self.axes[2].imshow(zoom)
        self.axes[2].figure.canvas.draw()

    def on_press_key(self, event):
        """Keyboards event.

        hit key list
        ------------
            v : -validate, close figs, and pass to next image,
            x : erase point mode. Hit x, if a point is in the vicinity it will be erase.
                Not reversable.
            d : enable/disable deplacement mode.
            'o' : original image for zoom axe
            'g' : gradient image for zoom axe
            'y' : gradient_y image for zoom axe
        """
        print(event.key)
        if event.key in ['v']:
            # pass to next image
            plt.close('all')

        elif event.key in ['x']:
            # erase closest point
            idx = 0
            while idx < len(self.xs):
                if (
                    abs(event.xdata - self.xs[idx]) < 8 and
                    abs(event.ydata - self.ys[idx]) < 8
                ):
                    break
                else:
                    idx += 1
            self.xs = np.delete(self.xs, idx)
            self.ys = np.delete(self.ys, idx)
            self.ys_interp = self.ys[np.argsort(self.xs)]
            self.xs_interp = np.sort(self.xs)
            self.interpolate()
            self.update_fig()

        elif event.key in ['d']:
            if self.dep_mode:  # si un point déjà selectionné
                if self.selected[0]:
                    self.dp1.set_data([[], []])
                    self.dp1.figure.canvas.draw()
                    self.selected[0] = False

                self.ax1.figure.canvas.mpl_disconnect(self.cid1)

                self.cid1 = self.line1.figure.canvas.mpl_connect(
                    'button_release_event', self.on_press
                )
                self.dep_mode = not self.dep_mode
                print('Add mode enabled' + '(Deplace mode disable)')
            else:
                self.ax1.figure.canvas.mpl_disconnect(self.cid1)
                self.cid1 = self.line1.figure.canvas.mpl_connect(
                    'button_release_event', self.deplace_release
                )
                self.dep_mode = not self.dep_mode
                print('Deplace mode enabled' + '(Add mode disable)')

        elif event.key in ['o', 'g', 'y']:
            if event.key == 'o':
                self.sel_ax_zoom = 0
            elif event.key == 'g':
                self.sel_ax_zoom = 1
            elif event.key == 'y':
                self.sel_ax_zoom = 3
            self.axes[2].figure.canvas.draw()

    def deplace_release(self, event):
        """
        Deplace a point.

        1 - Firts left clic gonna select a point.
        2 - Once a point is selected, other clic gonna deplace it.
        3 - Right clic to unselect the dot.
        4 - press 'd' to pass on add point mode.
        """
        if event.button == 3:
            """
            Clic droit pour déselectionner un point.
            Right clic to unselect a dot.
            """
            self.dp1.set_data([[], []])
            self.dp1.figure.canvas.draw()
            self.selected[0] = False
            return

        if event.button == 2:
            """
            Clic du milieu pour valider la figure.
            Middle clic to validate and pass through next figure.
            """
            plt.close()

        if self.selected[0] and event.button == 1:
            """There is a point selected."""
            x, y = self.dp1.get_xdata(), self.dp1.get_ydata()
            if not isinstance(x, (int, np.int, np.int64)):
                x = x[0]
                y = y[0]
            self.xs[self.selected[1]] = event.xdata
            self.ys[self.selected[1]] = event.ydata
            x = np.append(x, event.xdata)
            y = np.append(y, event.ydata)

            self.dp1.set_data(x, y)
            self.dp1.figure.canvas.draw()
            self.ys_interp = self.ys[np.argsort(self.xs)]
            self.xs_interp = np.sort(self.xs)
            self.interpolate()
            self.update_fig()

        elif not(self.selected[0]) and event.button == 1:
            """Select point to deplace."""
            idx = 0
            while idx < len(self.xs):
                if (
                    abs(event.xdata - self.xs[idx]) < 8 and
                    abs(event.ydata - self.ys[idx]) < 8
                ):
                    break
                else:
                    idx += 1
            if idx >= len(self.xs):
                print('Impossible to select a point in this area')
                self.selected[1] = []
            else:
                self.selected[1] = idx
                self.dp1.set_data(
                    self.xs[self.selected[1]],
                    self.ys[self.selected[1]]
                )
                self.dp1.figure.canvas.draw()
                self.selected[0] = True


class SplineBuilder_nonlinear(object):
    """
    Create Line on graph.

    Methods
    -------
    __call__
        on call the method plot the line

    """

    def __init__(self, ax1, ax2, ax3, ax4):
        """Do the initiation."""
        self.dep_mode = False  # at init no point is selected
        self.selected = [False, []]
        # first indicate if there is a selected point,
        # if True the second indicata is index

        self.values = ax2.get_images()[0].get_array()
        print(self.values)

        # line 1 #
        self.ax1 = ax1
        self.line1 = ax1.plot([], [])[0]
        self.line1.set_linestyle('-')
        self.line1.set_color('b')
        self.line1.set_linewidth(.8)

        self.pt1 = ax1.plot([], [])[0]
        print(self.pt1)
        self.pt1.set_linestyle('None')
        self.pt1.set_marker('o')
        self.pt1.set_markeredgecolor('b')

        self.dp1 = ax1.plot([], [])[0]
        print(self.dp1)
        self.pt1.set_linestyle('None')
        self.dp1.set_marker('s')
        self.dp1.set_markeredgecolor('r')
        self.dp1.set_markerfacecolor('none')

        # line 2 #
        self.line2 = ax2.plot([], [])[0]
        self.line2.set_linestyle('-')
        self.line2.set_color('b')
        self.line2.set_linewidth(.8)

        self.pt2 = ax2.plot([], [])[0]
        self.pt2.set_linestyle('None')
        self.pt2.set_marker('o')
        self.pt2.set_markeredgecolor('b')

        # line 3 #
        self.line3 = ax3.plot([], [])[0]
        self.line3.set_linestyle('-')
        self.line3.set_color('b')
        self.line3.set_linewidth(.8)

        self.pt3 = ax3.plot([], [])[0]
        self.pt3.set_linestyle('None')
        self.pt3.set_marker('o')
        self.pt3.set_markeredgecolor('b')

        # line 4 #
        self.line4 = ax4.plot([], [])[0]
        self.line4.set_linestyle('-')
        self.line4.set_color('b')
        self.line4.set_linewidth(.8)

        self.pt4 = ax4.plot([], [])[0]
        self.pt4.set_linestyle('None')
        self.pt4.set_marker('o')
        self.pt4.set_markeredgecolor('b')

        self.cid1 = self.ax1.figure.canvas.mpl_connect(
            'button_release_event', self.on_press
        )
        self.cid1_key = self.ax1.figure.canvas.mpl_connect(
            'key_release_event', self.on_press_key
        )

        self.xs = np.array(self.line1.get_xdata(), dtype=np.int)
        self.ys = np.array(self.line1.get_ydata(), dtype=np.int)
        self.xs_interp = np.array(self.line1.get_xdata(), dtype=np.int)
        self.ys_interp = np.array(self.line1.get_ydata(), dtype=np.int)

    def update_fig(self):
        """Update the figure."""
        self.pt1.set_data(self.xs, self.ys)
        self.pt2.set_data(self.xs, self.ys)
        self.pt3.set_data(self.xs, self.ys)
        self.pt4.set_data(self.xs, self.ys)

        self.line1.set_data(self.xs_interp, self.ys_interp)
        self.line2.set_data(self.xs_interp, self.ys_interp)
        self.line3.set_data(self.xs_interp, self.ys_interp)
        self.line4.set_data(self.xs_interp, self.ys_interp)

        self.line1.figure.canvas.draw()
        self.line2.figure.canvas.draw()
        self.line3.figure.canvas.draw()
        self.line4.figure.canvas.draw()

    def interpolate(self):
        """Construct the line interpolation."""
        if len(self.xs) > 2:

            rbf = interpolate.Rbf(self.xs_interp, self.ys_interp)
            self.xs_interp = np.arange(
                np.min(self.xs), np.max(self.xs), 1
            )
            self.ys_intepr = np.array(rbf(self.ys_interp), dtype=np.int)

            # tck = interpolate.splrep(self.xs_interp, self.ys_interp)
            # self.xs_interp = np.arange(
            #     np.min(self.xs), np.max(self.xs), 1
            # )
            # self.ys_interp = interpolate.splev(self.xs_interp, tck, der=0)

            # f = interpolate.interp1d(
            #     self.xs_interp, self.ys_interp,
            #     'cubic'
            # )
            # self.xs_interp = np.arange(
            #     np.min(self.xs), np.max(self.xs), 1
            # )
            # self.ys_interp = np.array(
            #     f(self.xs_interp), dtype=np.int
            # )
        elif len(self.xs) == 2:
            f = interpolate.interp1d(
                self.xs_interp, self.ys_interp,
                'slinear'
            )
            self.xs_interp = np.arange(
                np.min(self.xs), np.max(self.xs), 1
            )
            self.ys_interp = np.array(
                f(self.xs_interp), dtype=np.int
            )
    def on_press(self, event):
        """Create line on clicking on the graph.

        Actions
        -------
            Press left button to add a point,
            Middle button to finish the plot,
            Right to remove the last one.

        """
        zoom_cond = (
            self.line1.axes.get_navigate_mode() != 'ZOOM' or
            self.line2.axes.get_navigate_mode() != 'ZOOM' or
            self.line3.axes.get_navigate_mode() != 'ZOOM' or
            self.line4.axes.get_navigate_mode() != 'ZOOM'
        )
        if event.inaxes not in [
            self.line1.axes, self.line2.axes, self.line3.axes, self.line4.axes,
        ]:
            return

        elif event.button == 2:
            """
            Clic du milieu pour valider les points et
            passer à l'image suivante.
            """
            plt.close()

        elif (
            event.button == 3 and zoom_cond and
            len(self.xs) > 0 and len(self.ys) > 0
        ):
            """
            Clic droit supprime le dernier point.
            """
            if len(self.xs) > 0 and len(self.ys) > 0:
                self.xs = np.delete(self.xs, -1)
                self.ys = np.delete(self.ys, -1)

        elif event.button == 1 and zoom_cond:
            """
            Clic gauche, ajoute un point.
            Avec ce mode d'ajout, le code va retrier les points dans
            l'ordre des x croissant.
            Il ne sont donc pas prient en compte dans l'ordre d'apparition.
            """
            self.xs = np.append(self.xs, int(event.xdata))
            self.ys = np.append(self.ys, int(event.ydata))

            """
            Find the local maximum gradient

            ..todo

            """
            # vals_l = list(self.values[
            #     int(self.ys)-2:int(self.ys)+2,
            #     int(self.xs)-2:int(self.xs)+2
            # ])
            # line = vals_l.index(max(vals_l))
            # col = vals_l[line].index(max(vals_l[line]))
            # self.xs = self.xs + col-2
            # self.ys = self.ys + line-2

        self.ys_interp = self.ys[np.argsort(self.xs)]
        self.xs_interp = np.sort(self.xs)

        self.interpolate()
        self.update_fig()

    def on_press_key(self, event):
        """Keyboards event."""
        # print(event.key)
        if event.key in ['q', 'Q', 'escape', 'v', 'V', 'n', 'N', 'enter']:
            # pass to next image
            plt.close()
        elif event.key in ['x', 'X']:
            idx = 0
            while idx < len(self.xs):
                if (
                    abs(event.xdata - self.xs[idx]) < 8 and
                    abs(event.ydata - self.ys[idx]) < 8
                ):
                    break
                else:
                    idx += 1
            self.xs = np.delete(self.xs, idx)
            self.ys = np.delete(self.ys, idx)
            self.ys_interp = self.ys[np.argsort(self.xs)]
            self.xs_interp = np.sort(self.xs)
            self.interpolate()
            self.update_fig()

        elif event.key in ['c', 'C', 'd', 'D', 's', 'S']:
            if self.dep_mode:  # si un point déjà selectionné
                if self.selected[0]:
                    self.dp1.set_data([[], []])
                    self.dp1.figure.canvas.draw()
                    self.selected[0] = False

                self.ax1.figure.canvas.mpl_disconnect(self.cid1)

                self.cid1 = self.line1.figure.canvas.mpl_connect(
                    'button_release_event', self.on_press
                )
                self.dep_mode = not self.dep_mode
                print('Add mode enabled' + '(Deplace mode disable)')
            else:
                self.ax1.figure.canvas.mpl_disconnect(self.cid1)
                self.cid1 = self.line1.figure.canvas.mpl_connect(
                    'button_release_event', self.deplace_release
                )
                self.dep_mode = not self.dep_mode
                print('Deplace mode enabled' + '(Add mode disable)')

    def deplace_release(self, event):
        """
        Deplace a point.

        1 - Firts left clic gonna select a point.
        2 - Once a point is selected, other clic gonna deplace it.
        3 - Right clic to unselect the dot.
        4 - press 'd' to pass on add point mode.
        """
        if event.button == 3:
            """
            Clic droit pour déselectionner un point.
            Right clic to unselect a dot.
            """
            self.dp1.set_data([[], []])
            self.dp1.figure.canvas.draw()
            self.selected[0] = False
            return

        if event.button == 2:
            """
            Clic du milieu pour valider la figure.
            Middle clic to validate and pass through next figure.
            """
            plt.close()

        if self.selected[0] and event.button == 1:
            """There is a point selected."""
            x, y = self.dp1.get_xdata(), self.dp1.get_ydata()
            if not isinstance(x, (int, np.int, np.int64)):
                x = x[0]
                y = y[0]
            self.xs[self.selected[1]] = event.xdata
            self.ys[self.selected[1]] = event.ydata
            x = np.append(x, event.xdata)
            y = np.append(y, event.ydata)

            self.dp1.set_data(x, y)
            self.dp1.figure.canvas.draw()
            self.ys_interp = self.ys[np.argsort(self.xs)]
            self.xs_interp = np.sort(self.xs)
            self.interpolate()
            self.update_fig()

        elif not(self.selected[0]) and event.button == 1:
            """Select point to deplace."""
            idx = 0
            while idx < len(self.xs):
                if (
                    abs(event.xdata - self.xs[idx]) < 8 and
                    abs(event.ydata - self.ys[idx]) < 8
                ):
                    break
                else:
                    idx += 1
            if idx >= len(self.xs):
                print('Impossible to select a point in this area')
                self.selected[1] = []
            else:
                self.selected[1] = idx
                self.dp1.set_data(
                    self.xs[self.selected[1]],
                    self.ys[self.selected[1]]
                )
                self.dp1.figure.canvas.draw()
                self.selected[0] = True


class ContourLineBuilder(object):
    """
    Create Line on graph.

    Methods
    -------
    __call__
        on call the method plot the line

    """

    def __init__(self, ax1, ax2, geom):
        """Do the initiation."""
        self.dep_mode = False  # at init no point is selected
        self.selected = [False, []]
        # first indicate if there is a selected point,
        # if True the second indicata is index

        # get geometry values #
        self.rc_l = geom[0]
        self.rc_r = geom[1]
        self.rc = geom[2]
        self.zc = geom[3]
        self.z0 = geom[4]
        self.zf = geom[5]

        # Axes #
        self.ax1 = ax1  # Original image
        self.ax2 = ax2  # gradient image

        # lines #
        self.ax1_left_line = ax1.plot([], [])[0]
        self.ax1_left_line.set_linestyle('-')
        self.ax1_left_line.set_color('b')
        self.ax1_left_line.set_linewidth(.5)

        self.ax1_right_line = ax1.plot([], [])[0]
        self.ax1_right_line.set_linestyle('-')
        self.ax1_right_line.set_color('b')
        self.ax1_right_line.set_linewidth(.5)

        self.ax2_left_line = ax2.plot([], [])[0]
        self.ax2_left_line.set_linestyle('-')
        self.ax2_left_line.set_color('b')
        self.ax2_left_line.set_linewidth(.5)

        self.ax2_right_line = ax2.plot([], [])[0]
        self.ax2_right_line.set_linestyle('-')
        self.ax2_right_line.set_color('b')
        self.ax2_right_line.set_linewidth(.5)

        # points #
        self.ax1_left_pts = ax1.plot([], [])[0]
        self.ax1_left_pts.set_linestyle('None')
        self.ax1_left_pts.set_marker('+')
        self.ax1_left_pts.set_markeredgecolor('g')

        self.ax1_right_pts = ax1.plot([], [])[0]
        self.ax1_right_pts.set_linestyle('None')
        self.ax1_right_pts.set_marker('+')
        self.ax1_right_pts.set_markeredgecolor('r')

        self.ax1_top_pt = ax1.plot(
            self.rc, self.z0
        )[0]
        self.ax1_top_pt.set_linestyle('None')
        self.ax1_top_pt.set_marker('+')
        self.ax1_top_pt.set_markeredgecolor('r')
        self.ax1_top_pt.set_markerfacecolor('orange')

        self.ax2_left_pts = ax2.plot([], [])[0]
        self.ax2_left_pts.set_linestyle('None')
        self.ax2_left_pts.set_marker('+')
        self.ax2_left_pts.set_markeredgecolor('g')

        self.ax2_right_pts = ax2.plot([], [])[0]
        self.ax2_right_pts.set_linestyle('None')
        self.ax2_right_pts.set_marker('+')
        self.ax2_right_pts.set_markeredgecolor('r')

        self.ax2_top_pt = ax1.plot(
            [], []
        )[0]
        self.ax2_top_pt.set_linestyle('None')
        self.ax2_top_pt.set_marker('o')
        self.ax2_top_pt.set_markeredgecolor('r')
        self.ax2_top_pt.set_markerfacecolor('orange')

        # print ref points #
        self.ax1.plot(self.rc_l, self.zc, 'sr', mfc='none')
        self.ax2.plot(self.rc_l, self.zc, 'sr', mfc='none')

        self.ax1.plot(self.rc_r, self.zc, 'sr', mfc='none')
        self.ax2.plot(self.rc_r, self.zc, 'sr', mfc='none')

        self.ax1.plot(
            [self.rc, self.rc], [self.zc, self.zf],
            '--b', alpha=.5, linewidth=.5
        )
        self.ax2.plot(
            [self.rc, self.rc], [self.zc, self.zf],
            '--b', alpha=.5, linewidth=.5
        )

        # enable interactives modes #
        self.cid1 = self.ax1.figure.canvas.mpl_connect(
            'button_release_event', self.on_press
        )
        self.cid1_key = self.ax1.figure.canvas.mpl_connect(
            'key_release_event', self.on_press_key
        )

    def interpolate(self, xs, ys, side='left'):
        """Construct the line interpolation."""
        z_top = self.ax1_top_pt.get_ydata()
        if side == 'left':
            xs = np.insert(xs, 0, self.rc_l)
            ys = np.insert(ys, 0, self.zc)
            xs = np.append(xs, self.rc)
            ys = np.append(ys, z_top)

        elif side == 'right':
            xs = np.append(xs, self.rc_r)
            ys = np.append(ys, self.zc)
            xs = np.insert(xs, 0, self.rc)
            ys = np.insert(ys, 0, z_top)

        """
        # Polyfit method
        """
        p = np.polyfit(ys, xs, 3)  # x and y axes are intervated
        z = np.poly1d(p)
        y_int = np.array(np.arange(min(ys), max(ys), 1), dtype=np.int)
        x_int = z(y_int)

        """
        # UnivariateSpline method --> need to much data
        """
        # spl = UnivariateSpline(ys, xs)
        # # spl.set_smoothing_factor(0.5)
        # y_int = np.array(np.arange(min(ys), max(ys), 1), dtype=np.int)
        # x_int = spl(y_int)

        """
        # Spline cubic method.
        """
        xs = xs[np.argsort(ys)]
        ys = np.sort(ys)
        if len(xs) > 2:
            f = interpolate.interp1d(
                ys, xs,
                'quadratic'
            )
        elif len(xs) == 2:
            f = interpolate.interp1d(
                ys, xs,
                'slinear'
            )
        y_int = np.arange(
            np.min(ys), np.max(ys), 1
        )
        if len(xs) > 1:
            x_int = np.array(
                f(y_int), dtype=np.int
            )
        else:
            x_int = xs

        return x_int, y_int

    def on_press(self, event):
        """Create line on clicking on the graph.

        Actions
        -------
            Press left button to add a point,
            Middle button to finish the plot,
            Right to remove the last one.

        """
        zoom_cond = (
            self.ax1.get_navigate_mode() != 'ZOOM' or
            self.ax2.get_navigate_mode() != 'ZOOM'
        )
        if event.inaxes not in [
            self.ax1, self.ax2
        ]:
            return

        elif event.button == 2:
            """
            Clic du milieu pour valider les points et
            passer à l'image suivante.
            """
            plt.close()

        elif (
            event.button == 3 and zoom_cond and
            len(self.xs) > 0 and len(self.ys) > 0
        ):
            """
            Clic droit supprime le dernier point.
            """
            if len(self.xs) > 0 and len(self.ys) > 0:
                self.xs = np.delete(self.xs, -1)
                self.ys = np.delete(self.ys, -1)

        elif event.button == 1 and zoom_cond:
            """
            Clic gauche, ajoute un point.
            Avec ce mode d'ajout, le code va retrier les points dans
            l'ordre des x croissant.
            Il ne sont donc pas prient en compte dans l'ordre d'apparition.
            """
            if event.xdata < self.rc:
                xs = np.array(self.ax1_left_pts.get_xdata(), dtype=np.int)
                ys = np.array(self.ax1_left_pts.get_ydata(), dtype=np.int)
            else:
                xs = np.array(self.ax1_right_pts.get_xdata(), dtype=np.int)
                ys = np.array(self.ax1_right_pts.get_ydata(), dtype=np.int)

            if len(xs) < 7:
                xs = np.append(xs, int(event.xdata))
                ys = np.append(ys, int(event.ydata))

            if event.xdata < self.rc:
                self.ax1_left_pts.set_xdata(xs)
                self.ax1_left_pts.set_ydata(ys)

                self.ax2_left_pts.set_xdata(xs)
                self.ax2_left_pts.set_ydata(ys)

                x_int, y_int = self.interpolate(xs, ys, 'left')
                print('x_interpolation', x_int)
                print('y_interpolation', y_int)

                self.ax1_left_line.set_xdata(x_int)
                self.ax1_left_line.set_ydata(y_int)
                self.ax2_left_line.set_xdata(x_int)
                self.ax2_left_line.set_ydata(y_int)
            else:
                self.ax1_right_pts.set_xdata(xs)
                self.ax1_right_pts.set_ydata(ys)

                self.ax2_right_pts.set_xdata(xs)
                self.ax2_right_pts.set_ydata(ys)

                x_int, y_int = self.interpolate(xs, ys, 'right')

                self.ax1_right_line.set_xdata(x_int)
                self.ax1_right_line.set_ydata(y_int)
                self.ax2_right_line.set_xdata(x_int)
                self.ax2_right_line.set_ydata(y_int)

            self.ax1.figure.canvas.draw()
            # self.ax2.figure.canvas.draw()

    def on_press_key(self, event):
        """Keyboards event."""
        # print(event.key)
        if event.key in ['q', 'Q', 'escape', 'v', 'V', 'n', 'N', 'enter']:
            # pass to next image
            plt.close()
        elif event.key in ['x', 'X']:
            """To erase a point"""

            if event.xdata < self.rc:  # left side #
                xs = self.ax1_left_pts.get_xdata()
                ys = self.ax1_left_pts.get_ydata()
            else:  # right side #
                xs = self.ax1_right_pts.get_xdata()
                ys = self.ax1_right_pts.get_ydata()

            idx = 0
            while idx < len(xs):
                if (
                    abs(event.xdata - xs[idx]) < 8 and
                    abs(event.ydata - ys[idx]) < 8
                ):
                    break
                else:
                    idx += 1
            xs = np.delete(xs, idx)
            ys = np.delete(ys, idx)

            if event.xdata < self.rc:
                self.ax1_left_pts.set_xdata(xs)
                self.ax1_left_pts.set_ydata(ys)

                self.ax2_left_pts.set_xdata(xs)
                self.ax2_left_pts.set_ydata(ys)

                x_int, y_int = self.interpolate(xs, ys, 'left')

                self.ax1_left_line.set_xdata(x_int)
                self.ax1_left_line.set_ydata(y_int)

                self.ax2_left_line.set_xdata(x_int)
                self.ax2_left_line.set_ydata(y_int)
            else:
                self.ax1_right_pts.set_xdata(xs)
                self.ax1_right_pts.set_ydata(ys)

                self.ax2_right_pts.set_xdata(xs)
                self.ax2_right_pts.set_ydata(ys)

                x_int, y_int = self.interpolate(xs, ys, 'right')

                self.ax1_right_line.set_xdata(x_int)
                self.ax1_right_line.set_ydata(y_int)

                self.ax2_right_line.set_xdata(x_int)
                self.ax2_right_line.set_ydata(y_int)

            self.ax1.figure.canvas.draw()

        elif event.key in ['up', 'down']:
            """To change the position of top point"""
            z = int(self.ax1_top_pt.get_ydata())
            if event.key == 'up':
                z -= 1
                self.ax1_top_pt.set_ydata(z)
            else:
                z += 1
                self.ax1_top_pt.set_ydata(z)

            self.ax1.figure.canvas.draw()


    def deplace_release(self, event):
        """
        Deplace a point.

        1 - Firts left clic gonna select a point.
        2 - Once a point is selected, other clic gonna deplace it.
        3 - Right clic to unselect the dot.
        4 - press 'd' to pass on add point mode.
        """
        if event.button == 3:
            """
            Clic droit pour déselectionner un point.
            Right clic to unselect a dot.
            """
            self.dp1.set_data([[], []])
            self.dp1.figure.canvas.draw()
            self.selected[0] = False
            return

        if event.button == 2:
            """
            Clic du milieu pour valider la figure.
            Middle clic to validate and pass through next figure.
            """
            plt.close()

        if self.selected[0] and event.button == 1:
            """There is a point selected."""
            x, y = self.dp1.get_xdata(), self.dp1.get_ydata()
            if not isinstance(x, (int, np.int, np.int64)):
                x = x[0]
                y = y[0]
            self.xs[self.selected[1]] = event.xdata
            self.ys[self.selected[1]] = event.ydata
            x = np.append(x, event.xdata)
            y = np.append(y, event.ydata)

            self.dp1.set_data(x, y)
            self.dp1.figure.canvas.draw()
            self.ys_interp = self.ys[np.argsort(self.xs)]
            self.xs_interp = np.sort(self.xs)
            self.interpolate()
            self.update_fig()

        elif not(self.selected[0]) and event.button == 1:
            """Select point to deplace."""
            idx = 0
            while idx < len(self.xs):
                if (
                    abs(event.xdata - self.xs[idx]) < 8 and
                    abs(event.ydata - self.ys[idx]) < 8
                ):
                    break
                else:
                    idx += 1
            if idx >= len(self.xs):
                print('Impossible to select a point in this area')
                self.selected[1] = []
            else:
                self.selected[1] = idx
                self.dp1.set_data(
                    self.xs[self.selected[1]],
                    self.ys[self.selected[1]]
                )
                self.dp1.figure.canvas.draw()
                self.selected[0] = True
