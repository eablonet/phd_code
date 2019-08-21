# -*- coding: utf-8 -*-
"""
Create spline to track object.

Authors
-------
eablonet

"""


# import
# ------
import sys
import numpy as np
from scipy import interpolate

from scipy.interpolate import UnivariateSpline

from matplotlib import pyplot as plt


# usefull methods
# ---------------
def log(*x):
    """Print in sys."""
    print(*x, file=sys.stderr)


def euclidian_dist(x0, y0, x1, y1):
    """Euclidian distance."""
    d = np.sqrt((x0-x1)**2 + (y0-y1)**2)
    return d


class Point:
    """Specify the point attibutes."""

    def __init__(self):
        """Init the class."""
        self.x = []
        self.y = []

    def add_point(self, x, y):
        """Add a point by specifying location."""
        self.x.append(x)
        self.y.append(y)

    def remove_point(self, x, y):
        """Delete a point by specifying location."""
        for i, val in enumerate(self.x):
            if val == x and self.y[i] == y:
                del(self.x[i])
                del(self.y[i])

    def remove_last_point(self):
        """Remove last point in the list."""
        if len(self.x) > 0:
            del(self.x[-1])
            del(self.y[-1])

    def get_point(self):
        """Return points location."""
        return self.x, self.y

    def get_closer(self, x, y):
        """Return closer point."""
        d = 1920**4
        index = None
        for i in range(len(self.x)):
            d0 = euclidian_dist(self.x[i], self.y[i], x, y)
            if d0 < d:
                d = d0
                index = i
        return self.x[index], self.y[index]

    def get_sort_point(self):
        """Return sorted by x points locations."""
        y = [self.y[i] for i in np.argsort(self.x)]
        x = np.sort(self.x)
        return x, y

    def get_interpolation(self, kind='cubic'):
        """Return interpolate points."""
        xs, ys = self.get_sort_point()
        if len(self.x) in [0, 1]:
            return self.x, self.y
        elif len(self.x) in [2, 3]:
            x = np.arange(np.min(self.x), np.max(self.x), 1)
            f = interpolate.interp1d(xs, ys, kind='slinear')
            return x, f(x)
        else:
            x = np.arange(np.min(self.x), np.max(self.x), 1)
            f = interpolate.interp1d(xs, ys, kind=kind)
            return x, f(x)


class Image:
    """Special class to manipulate image."""
    def __init__(self):
        """Init the class."""
        self.im = None
        self.grad = None
        self.grady = None
        self.zoom = None

    def set_image(self, im):
        """Set the image."""
        self.im = im

    def set_grad(self, grad_im):
        """Set the gradient image."""
        self.grad = grad_im

    def set_grady(self, grady_im):
        """Set the image."""
        self.grady = grady_im

    def set_zoom(self, x, y, im='image', width=40, height=40):
        """Get zoom image.

        Is equivalent to crop image around a point.
        """
        case = {
            'image': self.im,
            'grad': self.grad,
            'grady': self.grady
        }
        w = int(width/2)
        h = int(height/2)
        im_w = np.shape(case[im])[1]
        im_h = np.shape(case[im])[0]

        x1 = (x-w)*(x > w)   # = x1-le si x1 > le 0 sinon 0
        y1 = (y-h)*(y > h)

        x2 = (
            (x+w) *
            (x < im_w - w+1) + im_w *
            (1 - (x < im_w - w+1))
        )
        y2 = (
            (y+h) *
            (y < im_h - h+1) + im_h *
            (1 - (y < im_h - h+1))
        )

        self.zoom = case[im][y1:y2, x1:x2]  # crop the image

        return x1, y1, x2, y2

    def get_image(self):
        """Return image."""
        return self.im

    def get_grad(self):
        """Return gradient image."""
        return self.grad

    def get_grady(self):
        """Return gradient_y image."""
        return self.grady

    def get_zoom(self):
        return self.zoom

class SplineBuilder(object):
    """
    Create Line on graph.

    Methods
    -------
    __call__
        on call the method plot the line

    """

    def __init__(self, image, image_grad, image_grady):
        """Do the initiation."""
        # store data
        # ----------
        self.image = Image()
        self.image.set_image(image)
        self.image.set_grad(image_grad)
        self.image.set_grady(image_grady)

        # Plot the main figures
        # --------------
        self.generate_home()

        # Allocate the point data
        # ----------------
        self.point = Point()

        # init states button
        # ------------------
        self.press_key = None
        self.zoom_center = []


        """Old stuff."""
        # self.intensity = intensity
        # self.z_intensity = z_intensity
        # self.r_intensity = r_intensity

        self.dep_mode = False  # at init no point is selected
        self.selected = [False, []]
        # first indicate if there is a selected point,
        # if True the second indicata is index

        # line 3 #
        # self.cross_x = ax3.plot(
        #     [], [], ls='-', color='tab:red', lw=2
        # )[0]  # x(cste)
        # self.cross_y = ax3.plot(
        #     [], [], ls='-', color='tab:red', lw=2
        # )[0]  # y(cste)

        # dots ???? kesako ?
        # -------------
        # self.dots = []
        # for a in self.axes2:
        #     self.dots.append(
        #         a.plot(
        #             [], [],
        #             ls='None',
        #             marker='o', ms=6, color='tab:red', mfc='None'
        #         )[0]
        #     )

        # enable connections
        # ------------------
        # self.cid1 = self.main_axes[0].figure.canvas.mpl_connect(
        #     'button_press_event', self.enable_press
        # )
        self.cid = self.main_fig.canvas.mpl_connect(
            'button_release_event', self.on_press
        )
        self.cid_key = self.main_fig.canvas.mpl_connect(
            'key_press_event', self.on_press_key
        )
        self.cid_release_key = self.main_fig.canvas.mpl_connect(
            'key_release_event', self.on_release_key
        )
        # self.cid1_move = self.main_axes[0].figure.canvas.mpl_connect(
        #     'motion_notify_event', self.on_motion
        # )

        # self.xs_interp = np.array(self.line1.get_xdata(), dtype=np.int)
        # self.ys_interp = np.array(self.line1.get_ydata(), dtype=np.int)

    def generate_home(self):
        """Generate the main figure."""
        self.main_fig = plt.figure(figsize=(16, 9))
        self.main_fig.canvas.set_window_title('Main Figure')
        self.main_fig.tight_layout(pad=.01, h_pad=.01, w_pad=.01)
        titles = [
            'Original image',
            'Gradient magnitude',
            'Zoom',
            'y-Gradient',
        ]

        self.main_axes = []
        self.main_frontLines = []
        self.main_frontPoints = []
        for i in range(4):
            self.main_axes.append(plt.subplot(2, 2, i+1))
            self.main_axes[i].axis('off')
            self.main_axes[i].set_title(titles[i])

            # front line
            self.main_frontLines.append(self.main_axes[i].plot([], [])[0])
            self.main_frontLines[i].set_linestyle('-')
            self.main_frontLines[i].set_color('b')
            self.main_frontLines[i].set_linewidth(.8)

            # clicked points
            self.main_frontPoints.append(self.main_axes[i].plot([], [])[0])
            self.main_frontPoints[i].set_linestyle('None')
            self.main_frontPoints[i].set_marker('o')
            self.main_frontPoints[i].set_markeredgecolor('b')

        self.main_axes[0].imshow(
            self.image.get_image(), cmap='pink', interpolation='bicubic'
        )
        self.main_axes[1].imshow(
            self.image.get_grad(), cmap='pink', interpolation='bicubic'
        )
        self.main_axes[3].imshow(
            self.image.get_grady(), cmap='pink', interpolation='bicubic'
        )

    def enable_press(self, event):
        self.press = True

    def update_fig(self):
        """Update the figure."""
        x, y = self.point.get_point()
        for pt in self.main_frontPoints:
            pt.set_data(x, y)

        x, y = self.point.get_interpolation()
        for line in self.main_frontLines:
            line.set_data(x, y)

        self.main_fig.canvas.draw()
        # self.line4.figure.canvas.draw()

    def _right_click(self):
        """Remove last Point."""
        self.point.remove_last_point()

    def _left_click(self, x, y):
        """Add a point."""
        self.point.add_point(x, y)

    def _middle_click(self, x, y):
        x, y = self.point.get_closer(x, y)
        self.point.remove_point(x, y)

    def on_press(self, event):
        """Create line on clicking on the graph.

        Actions
        -------
            Press left button to add a point,
            Middle button to remove closer point,
            Right to remove the last one.

        """
        # check if zoom is activated
        # --------------------------
        for line in self.main_frontLines:
            if line.axes.get_navigate_mode() == 'ZOOM':
                return

        # check if the click is not in one axe
        # ------------------
        if event.inaxes not in self.main_axes:
            return

        # action from click
        # -----------------
        if self.press_key is None:  # if we don't press a key, add/remove point
            x, y = int(event.xdata), int(event.ydata)
            if event.inaxes == self.main_axes[2]:
                x += self.zoom_center[0]
                y += self.zoom_center[1]

            if event.button == 1:
                # left click - add point
                self._left_click(x, y)

            elif event.button == 2:
                # middle click - remove closer
                self._middle_click(x, y)

            elif event.button == 3:
                # right click - remove last
                self._right_click()

        elif self.press_key == 'z':  # set zoom image
            axe = None
            for i in range(len(self.main_axes)):
                if event.inaxes == self.main_axes[i]:
                    axe = i
            if axe is None:  # cursor is no in one axe
                return

            case = {0: 'image', 1: 'grad', 3: 'grady'}
            x1, y1, _, _ = self.image.set_zoom(int(event.xdata), int(event.ydata), case[axe])
            self.zoom_center = [x1, y1]
            self.main_axes[2].imshow(
                self.image.get_zoom(),
                cmap='pink', interpolation='bicubic'
            )

        elif self.press_key == 'x':
            x, y = self.point.get_closer(int(event.xdata), int(event.ydata))
            self.point.remove_point(x, y)

        # update the figure
        # -----------------
        self.update_fig()

        # new quicker method
        """
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
        """

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
        axe = None
        for i in range(len(self.main_axes)):
            if event.inaxes == self.main_axes[i]:
                axe = i
        if axe is None:  # cursor is no in one axe
            return

        case = {0: 'image', 1: 'grad', 3: 'grady'}

        self.image.set_zoom(int(event.xdata), int(event.ydata), case[axe])

        self.main_axes[2].imshow(
            self.image.get_zoom(),
            cmap='pink', interpolation='bicubic'
        )

        self.main_fig.canvas.draw()

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
        if self.press_key is None:
            self.press_key = event.key


        """
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
    """

    def on_release_key(self, event):
        """When release a key."""
        if self.press_key == 'v':
            plt.close('all')

        self.press_key = None

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

    def add_line(self, x, y):
        """Add a ref line on the plot."""
        for a in self.main_axes:
            a.plot(
                x, y,
                ls='-.',
                color='r',
                alpha=.3
            )

    def add_points(self, x, y):
        """Add a ref line on the plot."""
        for a in self.main_axes:
            a.plot(
                x, y,
                ls='None', marker='x',
                mfc='None', mec='r', ms=4,
                alpha=.5,
            )


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

class SplineBuilderOld(object):
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
        plt.show()
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
