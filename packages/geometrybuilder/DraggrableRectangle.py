# -*- coding: utf-8 -*-
"""
Create a rectangle on an image.

Authors
-------
eablonet

"""
from matplotlib import pyplot as plt
import matplotlib.patches as patches


class DraggableRectangle(object):
    """Permit to drag a rectangle on a image."""

    def __init__(self, ax=None, rect=None):
        """Do the initiation."""
        if rect is None:
            self.rect = patches.Rectangle(
                (0, 0), 0, 0,
                linewidth=1, edgecolor='r', fill=None
            )
            ax.add_patch(self.rect)
        else:
            self.rect = rect

        self.xs = None
        self.ys = None
        self.xf = None
        self.yf = None
        self.press = None
        self.valid = False

        self.cidpress = self.rect.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        """Return positon value.

        Button Press
        ------------
        1 : left click
            store xs, ys (x_start, y_start)
        2 : middle click
            Valid the rectangle selection
        3 : right click
            no crop

        Parameters
        ----------
        event
            mouse event

        """
        if event.inaxes != self.rect.axes:
            return

        if event.button == 2:  # left click
            self.valid = True
            plt.close()

        elif event.button == 3:  # right click
            self.valid = True
            plt.close()

        elif event.button == 1:
            self.xs = int(event.xdata)
            self.ys = int(event.ydata)
            self.rect.set_x(self.xs)
            self.rect.set_y(self.ys)
            self.press = self.xs, self.ys

    def on_motion(self, event):
        """Display rectangle on motion mouse."""
        if self.press is None:
            return
        if event.inaxes != self.rect.axes:
            return

        self.rect.set_width(
            event.xdata - self.xs
        )
        self.rect.set_height(
            event.ydata - self.ys
        )

        self.rect.figure.canvas.draw()

    def on_release(self, event):
        """Return xf, yf on release."""
        if event.inaxes != self.rect.axes:
            self.xf = self.rect.get_height() + self.xs
            self.yf = self.rect.get_width() + self.ys
        else:
            self.xf = int(event.xdata)
            self.yf = int(event.ydata)

        self.rect.set_width(
            self.xf - self.xs
        )
        self.rect.set_height(
            self.yf - self.ys
        )

        self.press = None
        self.rect.figure.canvas.draw()

    def get_rect(self):
        """Return bottom left corner."""
        return min(self.xs, self.xf), max(self.xs, self.xf),\
            min(self.ys, self.yf), max(self.ys, self.yf)
