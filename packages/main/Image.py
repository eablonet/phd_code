# -*- coding: utf-8 -*-
"""
Image processing.

Allow you to define an object image with some useful properties.

:auth: Emeryk Ablonet - eablonet
:phD Student
:copyright: IMFT and Toulouse INP
:date : 2018

:version: 5.0
"""
import numpy as np
from scipy import signal

from matplotlib import pyplot as plt
from skimage.io import imread
from scipy.ndimage.interpolation import rotate

import clahe


class Image:
    """
    Class to treat an image.

    Can receive to type of file : '.tif' or '.npy'.

    .. todo :
        1. histImage show the histogram of the image
        2. Implement other type of file
    """

    def __init__(self, path=None):
        """
        Image Processing initiation.

        Method to define the imagepath.

        :param impath: is the path of the image to treat
        """
        self.path = None
        self.image = None
        self.width = None
        self.heigth = None
        self.size = None

        if path is not None:
            self.read_by_path(path)

    def __repr__(self):
        """Show info relative to the image."""
        txt = 'Image path : ' + self.path
        txt += '\nImage size (w x h) : ' + str(self.width)
        txt += 'x' + str(self.height)
        return txt

    def read_by_path(self, path):
        """Load the image by providing the path."""
        self.path = path
        if path[-4:] == '.tif':
            self.image = imread(path, as_gray=True)
        elif path[-4:] == '.npy':
            self.image = np.load(path)

        self.width = self.image.shape[1]
        self.height = self.image.shape[0]

        self.size = self.image.shape

    def set_image(self, image):
        """Load the image by providing the image itself."""
        self.path = None
        self.image = image
        self.width = self.image.shape[1]
        self.height = self.image.shape[0]

        self.size = self.image.shape

    def get_image(self):
        """Return the image."""
        return self.image

    def get_width(self):
        """Return the image width."""
        return self.width

    def get_height(self):
        """Return the image height."""
        return self.height

    def show_image(self, title=None):
        """Display the image."""
        plt.figure(figsize=(16, 9))
        if title is not None:
            plt.gcf().canvas.set_window_title(title)
        plt.imshow(self.image, cmap='gray')
        plt.show()

    def add_to_ax(self, ax, cmap='gray'):
        """Display the image on a given ax."""
        ax.imshow(self.image, cmap=cmap)

    def crop(self, xmin, xmax, ymin, ymax):
        """
        Permits to crop the image thanks to coordinate.

        input
        -----
            xmin, xmax: int
                Minimum and maximum positon to crop on x axis
            ymin, ymax: int
                Minimum and maximum positon to crop on y axis
        """
        self.image = self.get_crop(xmin, xmax, ymin, ymax)

        self.width = self.image.shape[1]
        self.height = self.image.shape[0]
        self.size = self.image.shape

    def get_crop(self, xmin, xmax, ymin, ymax):
        """Return image croped."""
        xmin = 0 if xmin < 0 else xmin
        ymin = 0 if ymin < 0 else ymin
        xmax = self.width if xmax > self.width else xmax
        ymax = self.height if ymax > self.height else ymax
        return self.image[ymin:ymax, xmin:xmax]

    def get_zoom(self, x, y, width=40, height=40):
        """Get zoom image.

        Is equivalent to crop image around a point.
        """
        w = int(width / 2)
        h = int(height/2)
        return self.get_crop(x-w, x+w, y-h, y+h)

    def set_zoom(self, x, y, width=40, height=40):
        """Set zoom image."""
        self.image = self.get_zoom(x, y, width, height)
        self.width = self.shape[1]
        self.height = self.shape[0]
        self.size = self.shape

    def get_col(self, row):
        """
        Return the column of a given row.

        input
        -----
            row: int
                Row of the colum to get
        """
        return self.image[:, row]

    def get_row(self, col):
        """
        Return the row of a given column.

        input
        -----
            col: int
                Column of the row to get
        """
        return self.image[col, :]

    def get_hist(self, bit=12):
        """Equalize automatic the histogram."""
        H, bins = np.histogram(self.image, 2**bit, [0, 2**bit])
        cs = H.cumsum()
        cs = cs * H.max() / cs.max()  # normalisation

        return H, cs

    def equalize_hist(self, bit=12):
        """Equalize automaticaly the image base on the histogram."""
        self.image = self.get_equalize(bit=bit)

    def get_equalize(self, bit=12):
        """Return the equalize automaticaly the image base on the histogram.

        The opencv equalizeHist is only working for 8bits. So we create our
        own equalizer using opencv theory :
        https://docs.opencv.org/2.4/modules/imgproc/doc/histograms.html?highlight=equalizehist#equalizehist
        1. calculate the nbit histogram H.
        2. Normalize it to nbit (sum of hist bins is 2**n)
        3. Compute the cumulative sum : cs(i) = ∑_j H(j) with j = [0:i]
        4. Transform the image using the new cumulative sum
        """
        H, bins = np.histogram(self.image, 2**bit, [0, 2**bit])

        cs = H.cumsum()
        cs = cs * H.max() / cs.max()  # normalisation
        cs = (cs - cs.min())*2**bit/(cs.max()-cs.min())  # linearisation
        return np.array(cs[self.image], dtype=np.uint16)

    def equalize_hist_by_clahe(self, limit=2, size=8):
        """Set the image equalize by clahe method.

        CLAHE : adaptative localized histogram equalizer.
        """
        self.image = self.get_equalize_by_clahe(limit=limit, size=size)

    def get_equalize_by_clahe(self, limit=2, size=8):
        """Return the equalize image by CLAHE (local enhancement).

        CLAHE : adaptative localized histogram equalizer.
        """
        return clahe.clahe(self.image, (size, size), limit)

    def set_gradient(self, size=3, type='sobel', out='mag', border='valid'):
        """Create the 2D derivative using."""
        self.image = self.get_gradient(
            size=3, type='sobel', out='mag', border='valid'
        )

    def get_gradient(self, size=3, type='sobel', out='mag', border='valid'):
        """Return gradient of the image."""
        if type == 'sobel':
            if size == 3:
                grady = np.array([
                    [-1, -2, -1],
                    [0, 0, 0],
                    [1, 2, 1]
                ]) / 8
                gradx = np.array([
                    [-1, 0, 1],
                    [-2, 0, 2],
                    [-1, 0, 1]
                ]) / 8
            elif size == 5:
                grady = np.array([
                    [-5, -8, -10, -8, -5],
                    [-4, -10, -20, -10, -4],
                    [0, 0, 0, 0, 0],
                    [4, 10, 20, 10, 4],
                    [5, 8, 10, 8, 5]
                ]) / 240
                gradx = np.array([
                    [-5, -4, 0, 4, 5],
                    [-8, -10, 0, 10, 8],
                    [-10, -20, 0, 20, 10],
                    [-8, -10, 0, 10, 8],
                    [-5, -4, 0, 4, 5]
                ]) / 240
        elif type == 'sobel-fedman':
            if size == 3:
                grady = np.array([
                    [-3, -10, -3],
                    [0, 0, 0],
                    [3, 10, 3]
                ]) / 32
                gradx = np.array([
                    [-3, 0, 3],
                    [-10, 0, 10],
                    [-3, 0, 3]
                ]) / 32
            elif size == 5:
                grady = np.array([
                    [-1, -2, -3, -2, -1],
                    [-1, -2, -6, -2, -1],
                    [0, 0, 0, 0, 0],
                    [1, 2, 6, 2, 1],
                    [1, 2, 3, 2, 1],
                ]) / 60
                gradx = np.array([
                    [-1, -1, 0, 1, 1],
                    [-2, -2, 0, 2, 2],
                    [-3, -6, 0, 6, 3],
                    [-2, -2, 0, 2, 2],
                    [-1, -1, 0, 1, 1]
                ]) / 60
        elif type == 'scharr':
            if size == 3:
                grady = np.array([
                    [-47, -162, -47],
                    [0, 0, 0],
                    [47, 162, 47]
                ])
                grady /= np.sum(np.abs(grady))
                gradx = np.array([
                    [-47, 0, 47],
                    [-162, 0, 162],
                    [-47, 0, 47]
                ])
                gradx /= np.sum(np.abs(gradx))

        gx = signal.convolve2d(self.image, gradx, border)
        gy = signal.convolve2d(self.image, grady, border)

        if out == 'x':
            return gx
        if out == 'y':
            return gy
        if out == 'mag':
            return np.sqrt(gx**2 + gy**2)
        if out == 'arg':
            return np.arctan2(gy, gx)

    def set_rotate(self, angle):
        """Rotate the image.

        input
        -----
            angle : float
                Angle of rotation in degree
        """
        self.image = self.get_rotate(angle)

    def get_rotate(self, angle):
        """Return rotated image.

        input
        -----
            angle : float
                Angle of rotation in degree
        """
        return rotate(self.image, angle, reshape=False)


if __name__ == '__main__':
    # create random image
    im = np.array(np.random.rand(200, 400)*255, dtype=np.int)
    im[100, 100:300] = 255
    im[50:150, 200] = 255
    im[90:110, 190:210] = 255

    # load image
    image = Image()
    image.set_image(im)
    # image.read_by_path('/path/to/image/')
    image.show_image('Original image')

    # localy equalize the image
    image.equalize_hist_by_clahe()
    image.show_image('Clahe image')

    # rotate the image
    image.set_rotate(30)
    image.show_image('30deg rotation')
    image.set_rotate(-30)
    image.show_image('-30deg rotation')

    # convolve the gradient
    image.set_gradient(size=3, type='scharr', out='mag')
    image.show_image('scharr gradient')
    image.set_gradient(size=3, type='sobel', out='mag')
    image.show_image('sobel gradient')
