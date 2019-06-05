# -*- coding: utf-8 -*-
"""
Create a rectangle on an image.

Authors
-------
eablonet

"""
from matplotlib.backends.qt_compat import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import (
    FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)
from matplotlib.figure import Figure


class ClaheWindow(QtWidgets.QMainWindow):
    """Choose Clahe parameters with UI."""

    def __init__(self, image):
        """Do the initiation.

        Parameter
        ---------
        image : StackProcessingObject
            provide the image to display in the app

        """
        super(ClaheWindow, self).__init__()
        self._main = QtWidgets.QWidget()
        self.image = image
        self.home()
        self.display()
        self.update()

    def home(self):
        """Do the home view."""
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)

        static_canvas = FigureCanvas(Figure(figsize=(12, 8)))
        self._static_ax = static_canvas.figure.subplots()

        layout.addWidget(static_canvas)
        self.addToolBar(NavigationToolbar(static_canvas, self))

        # slider limit #
        # ------------ #
        self.label_limit = QtWidgets.QLabel("Limit : 1.0")
        layout.addWidget(self.label_limit)

        h_layout = QtWidgets.QHBoxLayout()

        label_min = QtWidgets.QLabel("1")
        h_layout.addWidget(label_min)
        self.slider_limit = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider_limit.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.slider_limit.setTickInterval(1)
        self.slider_limit.setSingleStep(1)
        self.slider_limit.setMinimum(10)
        self.slider_limit.setMaximum(40)
        self.slider_limit.setValue(1)
        self.slider_limit.valueChanged.connect(self.value_changed)
        self.slider_limit.sliderReleased.connect(self.update)
        h_layout.addWidget(self.slider_limit)

        label_max = QtWidgets.QLabel("4")
        h_layout.addWidget(label_max)

        layout.addLayout(h_layout)

        # slider size #
        # ------------ #
        self.label_size = QtWidgets.QLabel("Size : (8x8)")
        layout.addWidget(self.label_size)

        h_layout = QtWidgets.QHBoxLayout()

        label_min = QtWidgets.QLabel("2")
        h_layout.addWidget(label_min)
        self.slider_size = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider_size.setTickPosition(QtWidgets.QSlider.TicksBothSides)
        self.slider_size.setTickInterval(4)
        self.slider_size.setSingleStep(1)
        self.slider_size.setMinimum(2)
        self.slider_size.setMaximum(32)
        self.slider_size.setValue(8)
        self.slider_size.valueChanged.connect(self.value_changed)
        self.slider_size.sliderReleased.connect(self.update)
        h_layout.addWidget(self.slider_size)

        label_max = QtWidgets.QLabel("32")
        h_layout.addWidget(label_max)

        layout.addLayout(h_layout)

    def display(self):
        """Display image."""
        self.imshow = self._static_ax.imshow(self.image.image, cmap='gray')

    def value_changed(self):
        """When value change."""
        self.label_limit.setText(
            "Limit : " + str(self.slider_limit.value()/10)
        )
        self.label_size.setText(
            "Size : (" +
            str(self.slider_size.value()) +
            "x" +
            str(self.slider_size.value()) +
            ")"
        )

    def update(self):
        """Update values."""
        size = self.slider_size.value()
        limit = self.slider_limit.value()/10
        im = self.image.equalize_hist_by_clahe(limit, size, output=True)
        self.imshow.set_data(im)
        self.imshow.set_clim(vmin=0, vmax=1)
        self.imshow.figure.canvas.draw()

    def get_data(self):
        """Return data.

        Returns
        -------
        limit : float
            Provide the limit to apply in clahe
        size : int
            Provide the size to apply in clahe

        """
        limit = self.slider_limit.value()/10
        size = self.slider_size.value()
        return limit, size
