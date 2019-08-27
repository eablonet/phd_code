# eablonet
# 2019
import numpy as np
import matplotlib.pyplot as plt


class Geometry(object):
    """
    Get rc_left, rc_right, zc, and z0/zf on figure.

    Methods
    -------
    __call__
        on call the method plot the line

    """

    def __init__(self, im_init, image_end, h, w):
        """Do the initiation."""
        self.ax1 = ax1
        self.h = h
        self.w = w

        # left line rc_left #
        self.rc_left_pos = int(w/4)
        self.rc_left = ax1.plot([], [], label='Left')[0]
        self.rc_left2 = ax2.plot([], [])[0]
        self.rc_left.set_linestyle('-')
        self.rc_left.set_color('g')
        self.rc_left.set_linewidth(.8)
        self.rc_left2.set_linestyle('-')
        self.rc_left2.set_color('g')
        self.rc_left2.set_linewidth(.8)

        # right line rc_right #
        self.rc_right_pos = int(3*w/4)
        self.rc_right = ax1.plot([], [], label='Right')[0]
        self.rc_right2 = ax2.plot([], [])[0]
        self.rc_right.set_linestyle('-')
        self.rc_right.set_color('r')
        self.rc_right.set_linewidth(.8)
        self.rc_right2.set_linestyle('-')
        self.rc_right2.set_color('r')
        self.rc_right2.set_linewidth(.8)

        # center line rc #
        self.rc_pos = int(w/2)
        self.rc = ax1.plot([], [], label='Center')[0]
        self.rc2 = ax2.plot([], [])[0]
        self.rc.set_linestyle('-')
        self.rc.set_color('b')
        self.rc.set_linewidth(.8)
        self.rc2.set_linestyle('-')
        self.rc2.set_color('b')
        self.rc2.set_linewidth(.8)

        # top line zf #
        self.zf_pos = int(h/4)
        self.zf = ax1.plot([], [], label='Top')[0]
        self.zf2 = ax2.plot([], [])[0]
        self.zf.set_linestyle('-')
        self.zf.set_color('c')
        self.zf.set_linewidth(.8)
        self.zf2.set_linestyle('-')
        self.zf2.set_color('c')
        self.zf2.set_linewidth(.8)

        # middle line z0 #
        self.z0_pos = int(h/2)
        self.z0 = ax1.plot([], [], label='Middle')[0]
        self.z02 = ax2.plot([], [])[0]
        self.z0.set_linestyle('-')
        self.z0.set_color('m')
        self.z0.set_linewidth(.8)
        self.z02.set_linestyle('-')
        self.z02.set_color('m')
        self.z02.set_linewidth(.8)

        # bottom line zc #
        self.zc_pos = int(3/4*h)
        self.zc = ax1.plot([], [], label='Bottom')[0]
        self.zc2 = ax2.plot([], [])[0]
        self.zc.set_linestyle('-')
        self.zc.set_color('y')
        self.zc.set_linewidth(.8)
        self.zc2.set_linestyle('-')
        self.zc2.set_color('y')
        self.zc2.set_linewidth(.8)

        self.ax1.legend(fancybox=True)

        self.update_fig()
        self.cid1 = self.ax1.figure.canvas.mpl_connect(
            'button_release_event', self.on_press
        )

        self.cid1_key = self.ax1.figure.canvas.mpl_connect(
            'key_release_event', self.on_press_key
        )

        self.sel = 'rc_left'

    def generate_home(self):
        self.main_fig = plt.figure()
        self.main_fig.canvas.set_window_title('Get the geometry')

        self.main_ax = []
        for i in range(2):
            self.main_ax.append(plt.subplot(1, 2, i+1))

    def update_fig(self):
        """Update the figure."""



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
        # print(event.key)
        switch = {
            'l': 'rc_left', 'r': 'rc_right', 'c': 'rc',
            't': 'zf', 'b': 'zc', 'm': 'z0'
        }

        if event.key in ['l', 'r', 't', 'b', 'c', 'm']:
            # changement de sélection
            self.sel = switch[event.key]

        elif event.key == 'up':  # cas up => z0 ou zc
            if self.sel == 'z0':
                self.z0_pos -= 1 if self.z0_pos > 0 else 0
                self.z0.set_data(
                    [0, self.w],
                    [self.z0_pos, self.z0_pos],
                )
                self.z02.set_data(
                    [0, self.w],
                    [self.z0_pos, self.z0_pos],
                )
            elif self.sel == 'zc':
                self.zc_pos -= 1 if self.zc_pos > 0 else 0
                self.zc.set_data(
                    [0, self.w],
                    [self.zc_pos, self.zc_pos],
                )
                self.zc2.set_data(
                    [0, self.w],
                    [self.zc_pos, self.zc_pos],
                )
            elif self.sel == 'zf':
                self.zf_pos -= 1 if self.zf_pos > 0 else 0
                self.zf.set_data(
                    [0, self.w],
                    [self.zf_pos, self.zf_pos],
                )
                self.zf2.set_data(
                    [0, self.w],
                    [self.zf_pos, self.zf_pos],
                )
        elif event.key == 'down':  # cas down => z0 ou zc
            if self.sel == 'z0':
                self.z0_pos += 1 if self.z0_pos < self.h else self.h
                self.z0.set_data(
                    [0, self.w],
                    [self.z0_pos, self.z0_pos],
                )
                self.z02.set_data(
                    [0, self.w],
                    [self.z0_pos, self.z0_pos],
                )
            elif self.sel == 'zc':
                self.zc_pos += 1 if self.zc_pos < self.h else self.h
                self.zc.set_data(
                    [0, self.w],
                    [self.zc_pos, self.zc_pos],
                )
                self.zc2.set_data(
                    [0, self.w],
                    [self.zc_pos, self.zc_pos],
                )
            elif self.sel == 'zf':
                self.zf_pos += 1 if self.zf_pos < self.h else self.h
                self.zf.set_data(
                    [0, self.w],
                    [self.zf_pos, self.zf_pos],
                )
                self.zf2.set_data(
                    [0, self.w],
                    [self.zf_pos, self.zf_pos],
                )

        elif event.key == 'left':  # cas left => rc_left ou rc_right
            if self.sel == 'rc_left':
                self.rc_left_pos -= 1 if self.rc_left_pos > 0 else 0
                self.rc_left.set_data(
                    [self.rc_left_pos, self.rc_left_pos],
                    [0, self.h],
                )
                self.rc_left2.set_data(
                    [self.rc_left_pos, self.rc_left_pos],
                    [0, self.h],
                )
            elif self.sel == 'rc_right':
                self.rc_right_pos -= 1 if self.rc_right_pos > 0 else 0
                self.rc_right.set_data(
                    [self.rc_right_pos, self.rc_right_pos],
                    [0, self.h],
                )
                self.rc_right2.set_data(
                    [self.rc_right_pos, self.rc_right_pos],
                    [0, self.h],
                )
            elif self.sel == 'rc':
                self.rc_pos -= 1 if self.rc_pos > 0 else 0
                self.rc.set_data(
                    [self.rc_pos, self.rc_pos],
                    [0, self.h],
                )
                self.rc2.set_data(
                    [self.rc_pos, self.rc_pos],
                    [0, self.h],
                )
        elif event.key == 'right':  # cas right => rc_left ou rc_right
            if self.sel == 'rc_left':
                self.rc_left_pos += 1 if self.rc_left_pos < self.w else self.w
                self.rc_left.set_data(
                    [self.rc_left_pos, self.rc_left_pos],
                    [0, self.h],
                )
                self.rc_left2.set_data(
                    [self.rc_left_pos, self.rc_left_pos],
                    [0, self.h],
                )
            elif self.sel == 'rc_right':
                self.rc_right_pos += 1 if self.rc_right_pos < self.w else self.w
                self.rc_right.set_data(
                    [self.rc_right_pos, self.rc_right_pos],
                    [0, self.h],
                )
                self.rc_right2.set_data(
                    [self.rc_right_pos, self.rc_right_pos],
                    [0, self.h],
                )
            elif self.sel == 'rc':
                self.rc_pos += 1 if self.rc_pos < self.w else self.w
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
