# -*- coding: utf-8 -*-
"""
Create spline to track object.

Authors
-------
eablonet

"""


# import
# ------
# sys
import sys

# math
import numpy as np
from scipy import interpolate
import peakutils as peaks  # find peaks method

# plot
from matplotlib import pyplot as plt
from matplotlib import patches

# local packages
# --------------
from packages.design import color as ea_color
from packages.library import main as ea_math


# usefull methods
# ---------------
def log(*x):
    """Print insys."""
    print(*x, file=sys.stderr)


def euclidian_dist(x0, y0, x1, y1):
    """Euclidian distance."""
    d = np.sqrt((x0 - x1)**2 + (y0 - y1)**2)
    return d


class Point:
    """Specify the point attibutes."""

    def __init__(self):
        """Init the class."""
        self.x = []
        self.y = []
        self.sx = []
        self.sy = []

    def add_point(self, x, y):
        """Add a point by specifying location.

        input
        -----
            x: int
                x position of the point
            y: int
                y position of the point
        """
        self.x.append(x)
        self.y.append(y)

    def set_new_location(self, x, y, x_new, y_new):
        """Redefine location of a point.

        input
        -----
            x, y: int or list
                location of the point to change
            x_new, y_new: int or list
                new coordinate of the point
        """
        # convert in list if int
        # --------------
        if type(x) is int:
            x, y, x_new, y_new = [x], [y], [x_new], [y_new]

        # loop to chane coordinate
        for i in range(len(x)):
            # find the index of x, y
            index = self.get_index_by_location(x[i], y[i])

            # change coordinate
            self.x[index] = x_new[i]
            self.y[index] = y_new[i]

    def set_selected_by_area(self, x0, y0, dx, dy):
        """Set point as seleted in an area.

        input
        -----
            xmin, ymin: int
                Top left coordinate of the area
            dx, dy: int
                width and height of the rectangle
        """
        # if dx/dy are negative
        # ------
        if dx < 0:
            x0 = x0 + dx
            dx *= -1
        if dy < 0:
            y0 = y0 + dy
            dy *= -1

        # loop over all ix to check if there are in the area
        # ------------
        for i in range(len(self.x)):
            if (
                self.x[i] < x0 + dx and self.x[i] > x0 and
                self.y[i] < y0 + dy and self.y[i] > y0
            ):
                # check if is already selected
                # -------------------
                isInSelection = False
                for j in range(len(self.sx)):
                    if self.sx[j] == self.x[i] and self.sy[j] == self.y[i]:
                        isInSelection = True

                # if not already selected -> add to selected list
                # ---------------
                if not isInSelection:
                    self.sx.append(self.x[i])
                    self.sy.append(self.y[i])

    def set_selected_by_location(self, x, y):
        """Set point as selected by location(s).

        input
        -----
            x, y: int or list
                location of the point(s) to select
        """
        if type(x) is int:  # convert in a list
            x = [x]
            y = [y]

        # loop for all given point
        # -----------
        for k in range(len(x)):

            # find the point in the points list
            # ----------------
            # if the point don't exist, the coordinate is ignore
            for i in range(len(self.x)):
                if self.x[i] == x[k] and self.y[i] == y[k]:
                    # check if is already selected
                    # -------------------
                    isInSelection = False
                    for j in range(len(self.sx)):
                        if self.sx[j] == x[k] and self.sy[j] == y[k]:
                            isInSelection = True

                    # if not already selected -> add to selected list
                    # ---------------
                    if not isInSelection:
                        self.sx.append(x[k])
                        self.sy.append(y[k])

    def remove_selected_by_area(self, x0, y0, dx, dy):
        """Remove point in an area.

        input
        -----
            xmin, ymin: int
                Top left coordinate of the area
            dx, dy: int
                width and height of the rectangle

        return
        ------
            c : int
                Number of deleted points
        """
        # if dx/dy are negative
        # --------------
        if dx < 0:
            x0 = x0 + dx
            dx *= -1
        if dy < 0:
            y0 = y0 + dy
            dy *= -1

        c = 0  # number of del point - avoid outof range
        for i in range(len(self.sx)):
            if (
                self.sx[i - c] < x0 + dx and self.sx[i - c] > x0 and
                self.sy[i - c] < y0 + dy and self.sy[i - c] > y0
            ):
                del(self.sx[i - c])
                del(self.sy[i - c])
                c += 1
        return c

    def reset_selection(self):
        """Reset selected point."""
        self.sx = []
        self.sy = []

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

    def get_one_point_by_location(self, x, y):
        """Return point of the coordinate.

        Return
        ------
            Point coordinate or empty if the point is not in the list.
        """
        for i in range(len(self.x)):
            if self.x[i] == x and self.y[i] == y:
                return self.x[i], self.y[i]

    def get_index_by_location(self, x, y):
        """Get the index for coordinate.

        return
        ------
            empty if the point does not exist. Else return the index.
        """
        for i in range(len(self.x)):
            if self.x[i] == x and self.y[i] == y:
                return i

    def get_point(self):
        """Return points location."""
        return self.x, self.y

    def get_selected_point(self):
        """Return selected point."""
        return self.sx, self.sy

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
        self.width = None
        self.height = None

    def set_image(self, im):
        """Set the image."""
        self.im = im
        self.width = np.shape(im)[1]
        self.height = np.shape(im)[0]
        print(self.width, self.height)

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
        w = int(width / 2)
        h = int(height / 2)
        im_w = np.shape(case[im])[1]
        im_h = np.shape(case[im])[0]

        x1 = (x - w) * (x > w)   # = x1-le si x1 > le 0 sinon 0
        y1 = (y - h) * (y > h)

        x2 = (
            (x + w) *
            (x < im_w - w + 1) + im_w *
            (1 - (x < im_w - w + 1))
        )
        y2 = (
            (y + h) *
            (y < im_h - h + 1) + im_h *
            (1 - (y < im_h - h + 1))
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
        """Return zoom image."""
        return self.zoom

    def get_height(self):
        """Return the height of the original image."""
        return self.height

    def get_width(self):
        """Return the width of the original image."""
        return self.width


class Rectangle:
    """Rectangle class to select."""

    def __init__(self):
        """Init the class."""
        self.x = 0
        self.y = 0
        self.dx = 0
        self.dy = 0
        self.rect = patches.Rectangle(
            (self.x, self.y), self.dx, self.dy,
            linewidth=2, edgecolor='tab:orange', fill=None,
            alpha=.5
        )

    def set_origin(self, x, y):
        """Set the origin point."""
        self.x = x
        self.y = y
        self.rect.set_x(x)
        self.rect.set_y(y)

    def set_size(self, x, y):
        """Set width and height of the rectangle."""
        self.dx = x - self.x
        self.dy = y - self.y

        self.rect.set_width(self.dx)
        self.rect.set_height(self.dy)

    def add_to_ax(self, ax):
        """Add the rectangle to an ax.

        caution
        -------
            The rectangle can be use only for one axe. To plot it in another
            please use reset before.
        """
        ax.add_patch(self.rect)

    def get_rect(self):
        """Return rectangle origin and dimension."""
        return self.x, self.y, self.dx, self.dy

    def reset(self):
        """Reset the rectangle."""
        self.rect.remove()
        self.__init__()


class Assistant:
    """Define the assistant class."""

    def __init__(self, image, r_guide):
        """Init the class."""
        self.r_guide = r_guide
        self.generate_home()

    def generate_home(self):
        """Create the main figure."""
        # init values
        # ------------
        self.intensity = []
        self.intensity_crop = []
        color = ea_color.color_met(7)

        self.main_fig = plt.figure(figsize=(16, 6))

        self.index = []
        for k in range(7):
            self.index.append([])
            self.intensity.append(
                ea_math.sliding_mean(
                    self.image.get_grad()[:, self.r_guide[k]], 3,
                )
            )

            front_r = int(self.ref_line[1][self.r_guide[k]])

            self.intensity_crop.append(
                self.intensity[k][front_r - 15:front_r + 6])

            self.main_axes.append(plt.subplot(1, 7, k + 1))
            self.main_axes[k].invert_yaxis()

            plt.plot(self.intensity_crop, range(21), '-', c=color[k])
            self.index[k].append(
                peaks.indexes(
                    self.intensity_crop,
                    thres=.3,
                    min_dist=3,
                )
            )

            plt.plot(
                [np.min(self.intensity_crop), np.max(self.intensity_crop)],
                [15, 15],
                '--', color='tab:red',
            )

            plt.plot(
                self.intensity_crop[self.index], self.index,
                ls='None',
                marker='o', ms=4,
                c='tab:green',
            )

    def get_point(self):
        """Return the index of point."""
        return self.index


class SplineBuilder(object):
    """
    Create Line on graph.

    Methods
    -------
    __call__
        on call the method plot the line

    """

    def __init__(self, image, image_grad, image_grady, geom):
        """Do the initiation."""
        self.help()

        # store data
        # ----------
        self.image = Image()
        self.image.set_image(image)
        self.image.set_grad(image_grad)
        self.image.set_grady(image_grady)
        self.geom = geom

        # Plot the main figures
        # --------------
        self.generate_home()

        # Allocate the point data
        # ----------------
        self.point = Point()

        # init states button
        # ------------------
        self.press_key = None
        self.rect = None
        self.drag = None
        self.guidLine = False
        self.ref_line = None

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
        self.p_click = self.main_fig.canvas.mpl_connect(
            'button_press_event', self.on_press_click
        )
        self.m_click = self.main_fig.canvas.mpl_connect(
            'motion_notify_event', self.on_motion_click
        )
        self.r_click = self.main_fig.canvas.mpl_connect(
            'button_release_event', self.on_release_click
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

    def help(self):
        """Print instruction."""
        print(
            """
            press key :
            -----------
                h : help
                s : select mode (enable/disable)
                g : guideline
                f : auto-find need to get
                d : deplacement mode (enable/disable)
                hold z + click : zoom on area
                hold x + click : remove closer point

            add mode
            --------
                left : add point
                right : remove last point
                middle : remove closer point

            select mode
            -----------
                s : disable select mode
                left-click + drag : create area, release to make the selection.
                right-click + drag : create area, release to remove point in
                    the area
                d : enable deplacement mode
                backspace : delete selected point

            deplace mode
            ------------
                d : disable deplace mode
                up/right/down/left : deplace all selected point (shift = x10)
                click-drage-release : deplace point
                s : enable select mode
            """
        )

    def generate_home(self):
        """Generate the main figure."""
        # Main figure cration
        # ---------------
        self.main_fig = plt.figure(figsize=(16, 9))
        self.main_fig.canvas.set_window_title('Main Figure')
        self.main_fig.tight_layout(pad=.01, h_pad=.01, w_pad=.01)

        # Disable shortcut
        # ----------------
        for param in plt.rcParams:
            if 'keymap' in param:
                plt.rcParams[param] = ''

        # Create the for axes
        # -------------------
        # Titles
        titles = [
            'Original image',
            'Gradient magnitude',
            'Zoom',
            'y-Gradient',
        ]

        # guideline
        color = ea_color.color_met(7)  # generate 7 colors
        self.r_guide = []  # init r_vector for guidLines \
        # for the 7 position where we gonna look at the intensity
        for i in range(7):
            self.r_guide.append(
                int(
                    i * 1.8 * self.geom.get_r0() / 6 +
                    self.geom.get_rc() - .9 * self.geom.get_r0()
                )
            )

        # init vectors objects
        # -------------------
        self.main_axes = []  # 4 main axes
        self.main_frontLines = []  # the front line for each axe
        self.main_frontPoints = []  # the clicked dot for each axe
        self.main_selectedPoints = []  # the selectedPoint for each axe
        self.main_guideLine = []  # the guidLine for each axe size = (4, 7)

        for i in range(4):
            self.main_axes.append(plt.subplot(2, 2, i + 1))
            self.main_axes[i].axis('off')
            self.main_axes[i].set_title(titles[i])

            # front line
            self.main_frontLines.append(self.main_axes[i].plot([], [])[0])
            self.main_frontLines[i].set_linestyle('-')
            self.main_frontLines[i].set_color('tab:blue')
            self.main_frontLines[i].set_linewidth(.8)

            # clicked points
            self.main_frontPoints.append(self.main_axes[i].plot([], [])[0])
            self.main_frontPoints[i].set_linestyle('None')
            self.main_frontPoints[i].set_marker('x')
            self.main_frontPoints[i].set_markeredgecolor('tab:blue')
            self.main_frontPoints[i].set_markerfacecolor('None')

            # selected points
            self.main_selectedPoints.append(self.main_axes[i].plot([], [])[0])
            self.main_selectedPoints[i].set_linestyle('None')
            self.main_selectedPoints[i].set_marker('x')
            self.main_selectedPoints[i].set_markeredgecolor('tab:orange')

            self.main_guideLine.append([])
            for k in range(7):
                self.main_guideLine[i].append(
                    self.main_axes[i].plot(
                        [], [], ls='--', c=color[k], alpha=.4, lw=2
                    )[0]
                )

        self.main_axes[0].imshow(
            self.image.get_image(), cmap='pink', interpolation='bicubic',
            # lut=2**12
        )
        self.main_axes[1].imshow(
            self.image.get_grad(), cmap='pink', interpolation='bicubic',
            # lut=2**12
        )
        self.main_axes[3].imshow(
            self.image.get_grady(), cmap='pink', interpolation='bicubic',
            # lut=2**12
        )

        self.disp_mode = self.main_axes[0].annotate(
            '',
            xy=(.975, .975), xycoords='figure fraction',
            horizontalalignment='right', verticalalignment='top',
            fontsize=12,
        )
        self.disp_mode.set_text('Mode : Add point')

    def update_fig(self):
        """Update the figure."""
        x, y = self.point.get_point()
        for pt in self.main_frontPoints:
            pt.set_data(x, y)

        x, y = self.point.get_interpolation()
        for line in self.main_frontLines:
            line.set_data(x, y)

        x, y = self.point.get_selected_point()
        for pt in self.main_selectedPoints:
            pt.set_data(x, y)
        self.main_fig.canvas.draw()

    def _right_click(self):
        """Remove last Point."""
        self.point.remove_last_point()

    def _left_click(self, x, y):
        """Add a point."""
        self.point.add_point(x, y)

    def _middle_click(self, x, y):
        x, y = self.point.get_closer(x, y)
        self.point.remove_point(x, y)

    def on_press_click(self, event):
        """Press click event."""
        # check if zoom is activated
        # --------------------------
        for line in self.main_frontLines:
            if line.axes.get_navigate_mode() == 'ZOOM':
                return

        # check if the click is not in one axe
        # ------------------
        if event.inaxes not in self.main_axes:
            return

        if self.press_key == 's':
            self.patch = self.rect.add_to_ax(event.inaxes)
            self.rect.set_origin(int(event.xdata), int(event.ydata))

    def on_motion_click(self, event):
        """When drag the mouse."""
        # check that there is no event in this mode
        # ----------------------------
        if self.press_key not in ['s']:
            return

        # check if zoom is activated
        # --------------------------
        for line in self.main_frontLines:
            if line.axes.get_navigate_mode() == 'ZOOM':
                return

        # check if the click is not in one axe
        # ------------------
        if event.inaxes not in self.main_axes:
            return

        if self.press_key == 's' and event.button is not None:
            self.rect.set_size(int(event.xdata), int(event.ydata))
            self.main_fig.canvas.draw()

    def on_release_click(self, event):
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

        # action from click/hold key
        # -------------------------
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
            x1, y1, _, _ = self.image.set_zoom(
                int(event.xdata),
                int(event.ydata),
                case[axe]
            )
            self.zoom_center = [x1, y1]
            self.main_axes[2].imshow(
                self.image.get_zoom(),
                cmap='pink', interpolation='bicubic'
            )

        elif self.press_key == 'x':  # remove closer point
            x, y = self.point.get_closer(int(event.xdata), int(event.ydata))
            self.point.remove_point(x, y)

        elif self.press_key == 's':  # select points
            self.rect.set_size(int(event.xdata), int(event.ydata))
            if event.button == 1:
                self.point.set_selected_by_area(*self.rect.get_rect())
            elif event.button == 3:
                self.point.remove_selected_by_area(*self.rect.get_rect())
            self.rect.reset()
            self.update_fig()

        elif self.press_key == 'd':  # deplace point
            None

        # update the figure
        # -----------------
        self.update_fig()

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
            for z in self.z_intensity[y_curs - 15:y_curs + 15]:
                dx = x_curs - self.intensity[z, n]
                dy = y_curs - z
                if np.sqrt(dx**2 + dy**2) < d:
                    d = np.sqrt(dx**2 + dy**2)
                    z_out = z

            # find the local maximum [+- 4px]
            it_max = 0  # init intensity maximum
            z_max = 0  # init z of maximum intensity
            for z in self.z_intensity[z_out - 4:z_out + 4]:
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

        elif self.press_key == 's':
            if self.rect is None:  # enable select mode
                self.rect = Rectangle()
                self.disp_mode.set_text(
                    """Mode : Select/remove point (left/right click)"""
                )
                self.update_fig()

            elif event.key in ['s', 'd']:  # disable select mode / change mode
                self.press_key = None
                self.rect = None
                self.disp_mode.set_text('Mode : Add point')
                if event.key == 'd':
                    self.press_key = event.key
                    self.drag = [0, 0]
                    self.disp_mode.set_text(
                        """Mode : Deplace point (arrow or drag-click)"""
                    )
                self.update_fig()

            elif event.key == 'backspace':  # erase selected point
                x, y = self.point.get_selected_point()
                for i in range(len(x)):
                    self.point.remove_point(x[i], y[i])
                    self.point.reset_selection()
                    self.update_fig()

        elif self.press_key == 'd':
            if self.drag is None:  # enable deplacement mode
                self.drag = [0, 0]
                self.disp_mode.set_text(
                    """Mode : Deplace point (arrow or drag-click)"""
                )
                self.update_fig()

            elif event.key in ['d', 's']:  # disable deplacement mode / change
                # mode
                self.drag = None
                self.press_key = None
                self.disp_mode.set_text('Mode : Add point')
                if event.key == 's':
                    self.press_key = event.key
                    self.rect = Rectangle()
                    self.disp_mode.set_text(
                        """Mode : Select/remove point (left/right click)"""
                    )
                self.update_fig()

            elif len(self.point.get_selected_point()) < 1:
                print(
                    """To deplace point, please select one or
                     more point before."""
                )

            elif 'up' in event.key:
                x, y = self.point.get_selected_point()
                ynew = [pt - (10 if 'shift' in event.key else 1) for pt in y]
                self.point.set_new_location(x, y, x, ynew)
                self.point.reset_selection()
                self.point.set_selected_by_location(x, ynew)
                self.update_fig()

            elif 'down' in event.key:
                x, y = self.point.get_selected_point()
                ynew = [pt + (10 if 'shift' in event.key else 1) for pt in y]
                self.point.set_new_location(x, y, x, ynew)
                self.point.reset_selection()
                self.point.set_selected_by_location(x, ynew)
                self.update_fig()

            elif 'right' in event.key:
                x, y = self.point.get_selected_point()
                xnew = [pt + (10 if 'shift' in event.key else 1) for pt in x]
                self.point.set_new_location(x, y, xnew, y)
                self.point.reset_selection()
                self.point.set_selected_by_location(xnew, y)
                self.update_fig()

            elif 'left' in event.key:
                x, y = self.point.get_selected_point()
                xnew = [pt - (10 if 'shift' in event.key else 1) for pt in x]
                self.point.set_new_location(x, y, xnew, y)
                self.point.reset_selection()
                self.point.set_selected_by_location(xnew, y)
                self.update_fig()

        elif self.press_key == 'h':
            self.help()
            self.press_key = None

        elif self.press_key == 'g':
            if not self.guidLine:
                for i in range(4):
                    for k in range(7):
                        self.main_guideLine[i][k].set_data(
                            [self.r_guide[k], self.r_guide[k]],
                            [0, self.image.get_height()]
                        )

            else:
                for i in range(4):
                    for k in range(7):
                        self.main_guideLine[i][k].set_data([], [])

            self.update_fig()
            self.guidLine = not self.guidLine
            self.press_key = None

        elif self.press_key == 'f':
            self._assistant()
            self.press_key = None

        else:
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

    def _assistant(self):
        intensity_r = []
        plt.figure(figsize=(16, 6))

        color = ea_color.color_met(7)

        for k in range(7):
            intensity_r.append(
                ea_math.sliding_mean(
                    self.image.get_grad()[:, self.r_guide[k]], 3,
                )
            )

            front_r = int(self.ref_line[1][self.r_guide[k]])

            intensity_r_crop = intensity_r[k][front_r - 15:front_r + 6]

            ax = plt.subplot(1, 7, k + 1)
            ax.invert_yaxis()

            plt.plot(intensity_r_crop, range(21), '-', c=color[k])
            index = peaks.indexes(
                intensity_r_crop,
                thres=.3,
                min_dist=3,
            )

            plt.plot(
                [np.min(intensity_r_crop), np.max(intensity_r_crop)],
                [15, 15],
                '--', color='tab:red',
            )

            plt.plot(
                intensity_r_crop[index], index,
                ls='None',
                marker='o', ms=4,
                c='tab:green',
            )
            plt.show()

    def add_line(self, x, y):
        """Add a ref line on the plot."""
        for a in self.main_axes:
            a.plot(
                x, y,
                ls='-.',
                color='r',
                alpha=.3
            )

    def set_ref_line(self, x, y):
        """Add a ref line on the plot."""
        # extend the line on the entire width of the figure
        # ------
        # we assume that x is sorted
        xnew = np.arange(self.image.get_width())
        ynew = np.zeros(self.image.get_width())
        if len(x) > 0:
            ynew[x] = y
            ynew[xnew < x[0]] = y[0]
            ynew[xnew > x[-1]] = y[-1]

        self.ref_line = [xnew, ynew]
        for a in self.main_axes:
            a.plot(
                xnew, ynew,
                ls='-.',
                color='tab:green',
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


class Geometry(object):
    """
    Get rc_left, rc_right, zc, and z0/zf on figure.

    Methods
    -------
    __call__
        on call the method plot the line

    """

    def __init__(self, im_init, im_end):
        """Do the initiation."""
        # import data
        # ----------
        self.image = Image()
        self.image.set_image(im_init)
        self.image.set_grad(im_end)

        # init data
        # ---------
        self.rl = None
        self.rc = None
        self.rr = None
        self.zf = None
        self.z0 = None
        self.zb = None

        # init button and selection
        # ----------------
        self.n_line_selected = 0

        # generate main figure
        # --------------------
        self.generate_home()

        # connect button/click
        # --------------------
        # self.main_fig.canvas.mpl_connect(
        #     'button_release_event', self.on_press
        # )

        self.main_fig.canvas.mpl_connect(
            'key_release_event', self.on_press_key
        )

    def generate_home(self):
        """Generate the main ux."""
        # create figure
        # -----------
        self.main_fig = plt.figure(figsize=(16, 9))
        self.main_fig.canvas.set_window_title('Get the geometry')

        # get width and height
        # -----------
        w, h = self.image.get_width(), self.image.get_height()

        # generate color and vectors position
        # ------------------
        color = ea_color.color_met(5)
        self.x = [
            [int(w/4), int(w/4)],
            [int(w/2), int(w/2)],
            [int(3*w/4), int(w/4)],
            *[[0, w]]*3,
        ]
        self.y = [
            *[[0, h]]*3,
            [int(h/4), int(h/4)],
            [int(h/2), int(h/2)],
            [int(3*h/4), int(h/4)],
        ]

        # loop to create axes, lines and annotation
        self.main_ax = []
        self.main_lines = []
        self.main_textLine = []
        for i in range(2):
            self.main_ax.append(plt.subplot(1, 2, i+1))
            self.main_ax[i].axis('off')

            # loop over lines
            self.main_lines.append([])
            self.main_textLine.append([])
            for k in range(6):
                self.main_lines.append(
                    self.main_ax[i].plot(
                        self.x[k], self.y[k],
                        ls='-', lw='.8',
                        c=color[k],
                    )[0]
                )
                self.main_textLine.append(
                    self.main_ax[i].annotate(
                        str(k),
                        xy=(self.x[k][0], self.y[k][0])
                    )
                )

        self.main_ax[0].imshow(
            self.image.get_image(),
            cmap='pink'
        )
        self.main_ax[1].imshow(
            self.image.get_grad(),
            cmap='pink'
        )

    def update_fig(self):
        """Update the figure."""
        for i in range(2):
            for k in range(6):
                self.main_lines[i][k].set_data(
                    self.x[k], self.y[k]
                )
                self.main_textLine.set_xy(
                    self.x[k][0], self.y[k][0]
                )
        self.main_fig.canvas.draw()

    def on_press(self, event):
        """Create line on clicking on the graph.

        Actions
        -------
            Press left button to add a point,
            Middle button to finish the plot,
            Right to remove the last one.

        """
        zoom_cond = (
            self.zc.axes.get_navigate_mode() != 'ZOOM',
            self.zc2.axes.get_navigate_mode() != 'ZOOM'
        )
        if event.inaxes not in [self.zc.axes, self.zc2.axes]:
            return

        elif event.button == 2:
            """
            Clic du milieu pour valider les points et
            passer à l'image suivante.
            """
            plt.close()

        elif event.button == 1 and zoom_cond:
            """
            Clic gauche, ajoute un point.
            Avec ce mode d'ajout, le code va retrier les points dans
            l'ordre des x croissant.
            Il ne sont donc pas prient en compte dans l'ordre d'apparition.
            """
            x = int(event.xdata)
            y = int(event.ydata)

            if self.sel == 'z0':
                self.z0_pos = y
                self.z0.set_data(
                    [0, self.w],
                    [self.z0_pos, self.z0_pos],
                )
                self.z02.set_data(
                    [0, self.w],
                    [self.z0_pos, self.z0_pos],
                )
            elif self.sel == 'zc':
                self.zc_pos = y
                self.zc.set_data(
                    [0, self.w],
                    [self.zc_pos, self.zc_pos],
                )
                self.zc2.set_data(
                    [0, self.w],
                    [self.zc_pos, self.zc_pos],
                )
            elif self.sel == 'zf':
                self.zf_pos = y
                self.zf.set_data(
                    [0, self.w],
                    [self.zf_pos, self.zf_pos],
                )
                self.zf2.set_data(
                    [0, self.w],
                    [self.zf_pos, self.zf_pos],
                )
            elif self.sel == 'rc_left':
                self.rc_left_pos = x
                self.rc_left.set_data(
                    [self.rc_left_pos, self.rc_left_pos],
                    [0, self.h],
                )
                self.rc_left2.set_data(
                    [self.rc_left_pos, self.rc_left_pos],
                    [0, self.h],
                )
            elif self.sel == 'rc_right':
                self.rc_right_pos = x
                self.rc_right.set_data(
                    [self.rc_right_pos, self.rc_right_pos],
                    [0, self.h],
                )
                self.rc_right2.set_data(
                    [self.rc_right_pos, self.rc_right_pos],
                    [0, self.h],
                )
            elif self.sel == 'rc':
                self.rc_pos = x
                self.rc.set_data(
                    [self.rc_pos, self.rc_pos],
                    [0, self.h],
                )
                self.rc2.set_data(
                    [self.rc_pos, self.rc_pos],
                    [0, self.h],
                )

        # re-draw figure
        self.zc.figure.canvas.draw()

    def on_press_key(self, event):
        """Keyboards event."""
        if event.key in ['0', '1', '2', '3', '4', '5']:
            # changement de sélection among 0, 1, 2, 3, 4, 5
            self.n_line_selected = int(event.key)

        new_pos = 10 if 'shift' in event.key else 1
        if 'up' in event.key or 'left' in event.key:
            new_pos *= -1

        cond = False
        for w in event.key.split('-'):
            if w in ['up', 'down', 'left', 'right']:
                cond = True

        if cond:
            if self.n_line_selected < 3:
                self.x[self.n_line_selected] = [
                    self.x[self.n_line_selected][0] + new_pos,
                    self.x[self.n_line_selected][0] + new_pos
                ]
            else:
                self.y[self.n_line_selected] = [
                    self.y[self.n_line_selected][0] + new_pos,
                    self.y[self.n_line_selected][0] + new_pos
                ]

            self.update_fig()

    def get_geom(self):
        """Return dictionary containing geometrical data."""
        geom = {
            'rl': self.rl,
            'rr': self.rr,
            'rc': self.rc,
            'z0': self.z0,
            'zf': self.zf,
            'zb': self.zb,
        }
        return geom

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
                x_min = x_sort[(x_sort == x_new) - 1]
                x_max = x_sort[(x_sort == x_new) + 1]

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
            for z in self.z_intensity[y_curs - 15:y_curs + 15]:
                dx = x_curs - self.intensity[z, n]
                dy = y_curs - z
                if np.sqrt(dx**2 + dy**2) < d:
                    d = np.sqrt(dx**2 + dy**2)
                    z_out = z

            # find the local maximum [+- 4px]
            it_max = 0  # init intensity maximum
            z_max = 0  # init z of maximum intensity
            for z in self.z_intensity[z_out - 4:z_out + 4]:
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

        x1 = (x1 - le) * (x1 > le)   # = x1-le si x1 > le 0 sinon 0
        y1 = (y1 - le) * (y1 > le)

        x2 = (x2 + le) * (x2 < np.shape(self.image)[1] - le + 1) + np.shape(
            self.image)[1] * (1 - (x2 < np.shape(self.image)[1] - le + 1))
        y2 = (y2 + le) * (y2 < np.shape(self.image)[0] - le + 1) + np.shape(
            self.image)[0] * (1 - (y2 < np.shape(self.image)[0] - le + 1))

        zoom = (
            self.image * (self.sel_ax_zoom == 0) +
            self.grad * (self.sel_ax_zoom == 1) +
            self.grady * (self.sel_ax_zoom == 3)
        )  # select which image to use (change pressing 'o, 'g or 'y')

        zoom = zoom[y1:y2, x1:x2]  # crop the image

        xc = int(np.shape(zoom)[1] / 2)  # x_center_zoom
        yc = int(np.shape(zoom)[0] / 2)  # y_center_zoom

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
