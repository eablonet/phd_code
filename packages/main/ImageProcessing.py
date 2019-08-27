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

from skimage.filters import threshold_otsu, threshold_local, gaussian
from skimage.morphology import reconstruction
from skimage.exposure import equalize_adapthist

import cv2 as cv


class ImageProcessing():
    """
    Class to treat an image.

    Can receive to type of file : '.tif' or '.npy'.

    .. todo :
        1. histImage show the histogram of the image
        2. Implement other type of file
    """

    def __init__(self, impath):
        """
        Image Processing initiation.

        Method to define the imagepath.

        :param impath: is the path of the image to treat
        """
        self.image_path = impath
        if impath[-4:] == '.tif':
            self.image = imread(impath, as_gray=True)
        elif impath[-4:] == '.npy':
            self.image = np.load(impath)

        self.width = self.image.shape[1]
        self.height = self.image.shape[0]

        self.size = self.image.shape

    def __repr__(self):
        """
        Show info relative to the image.

        Print :
            1. directory which is loaded
            2. number of image loaded
            3. current image selected
            4. first image number of the stack
            5. last image number of the stack
        """
        print('Image path : ', self.image_path)
        print('Image size (w x h) : ', self.widths + 'x' + self.height)

    def get_image(self):
        """Return the image."""
        return self.image

    def get_width(self):
        """Return the image width."""
        return self.width

    def get_height(self):
        """Return the image height."""
        return self.height

    def show_image(self):
        """Display the image."""
        plt.figure(figsize=(16, 9))
        plt.imshow(self.image, cmap='gray')
        plt.show()

    def crop(self, xmin, xmax, ymin, ymax):
        """
        Permits to crop the image thanks to coordinate.

        :param xmin: mininaml x coordinate
        :param xmax: maximal x coordinate
        :param ymin: minimal y coordinate
        :param ymax: maximal y coordinate
        """
        if xmin < 0:
            print('xmin might be at least 0 - 0 will be used')
            xmin = 0
        if ymin < 0:
            print('ymin might be at least 0 - 0 will be used')
            ymin = 0
        if xmax > self.size[1]:
            print(
                'xmax might be less than', str(self.size[1]),
                '- this value will be used'
            )
            xmax = self.size[1]
        if ymax > self.size[0]:
            print(
                'ymax might be less than', str(self.size[0]),
                '- this value will be used'
            )
            ymax = self.size[0]
        self.image = self.image[
                ymin:ymax,
                xmin:xmax
            ]
        self.size = self.image.shape

    def get_crop(self, xmin, xmax, ymin, ymax):
        """
        Permits to crop the image thanks to coordinate.

        input:
            xmin: int
                mininaml x coordinate
            xmax: int
                maximal x coordinate
            ymin: int
                minimal y coordinate
            ymax: int
                maximal y coordinate
        """
        if xmin < 0:
            xmin = 0
        if ymin < 0:
            ymin = 0
        if xmax > self.size[1]:
            xmax = self.size[1]
        if ymax > self.size[0]:
            ymax = self.size[0]

        image = self.image[
                ymin:ymax,
                xmin:xmax
            ]

        return image

    def row_intensity(self, row, N=7):
        """
        Display the intensity on designed lined number.

        :param row: row number
        :param N: length for average moving calcul
        """
        intensity = self.image[row, :]
        it_a = np.convolve(
            intensity, np.ones((N,))/N, mode='valid'
        )
        return intensity, it_a

    def col_intensity(self, col, N=7):
        """
        Display the intensity on designed lined number.

        :param row: row number
        :param N: length for average moving calcul
        """
        intensity = self.image[:, col]
        it_a = np.convolve(intensity, np.ones((N,))/N, mode='valid')
        return intensity, it_a

    def plot_hist(self, bit=12):
        """Equalize automatic the histogram."""
        H, bins = np.histogram(self.image, 2**bit, [0, 2**bit])
        cs = H.cumsum()
        cs = cs * H.max() / cs.max()  # normalisation
        """
        plt.figure()
        plt.plot(H, color='r')
        plt.plot(cs, color='b')
        plt.show()
        """
        return H, cs

    def equalize_hist(self, bit=12):
        """Equalize automatic the histogram.

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
        self.image = np.array(cs[self.image], dtype=np.uint16)

    def equalize_hist_by_clahe(self, limit=2, size=8, output=False):
        """Equalize automatic the histogram by local method.

        CLAHE : adaptative localized histogram equalizer.
        """
        # opencv method
        # clahe = cv.createCLAHE(clipLimit=limit, tileGridSize=(size, size))
        # self.image = clahe.apply(self.image)

        # skiimage method
        """
        self.image = np.array(
            equalize_adapthist(
                self.image, clip_limit=limit, nbins=2**12
            )*2**12, dtype=np.int
        )
        """

        if not output:
            self.image = clahe.clahe(self.image, (size, size), limit)
        else:
            return clahe.clahe(self.image, (size, size), limit)

    def gradient(self, size=3, type='sobel', out='mag', border='valid', im=None):
        """Create the 2D derivative using ..type method."""
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
        elif type == 'scharr':
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

        if im is None:
            gy = signal.convolve2d(self.image, grady, border)
            gx = signal.convolve2d(self.image, gradx, border)
        else:
            gy = signal.convolve2d(im, grady, border)
            gx = signal.convolve2d(im, gradx, border)

        if out == 'x':
            return gx
        if out == 'y':
            return gy
        if out == 'mag':
            gmag = np.sqrt(gy**2 + gx**2)
            return gmag
        if out == 'arg':
            garg = np.arctan2(gy, gx)
            return garg

    def rotate(self, angle, im=None):
        """Rotate the image.

        Parameters
        ----------
        angle : float
            Angle of rotation in degree

        """
        if im is None:
            rot = rotate(self.image, angle, reshape=False)
        else:
            rot = rotate(im, angle, reshape=False)

        return rot


if __name__ == '__main__':
    from skimage.filters.rank import entropy
    from skimage.morphology import disk
    path = \
        '/Users/eablonet/Documents/0_phd/2_results/export_image/rep1/rep1n06/\
lastOp/i_0010.npy'
    current_image = ImageProcessing(path)
    current_image.equalize_hist_by_clahe(.2, 8)
    # current_image.gradient(type='scharr', size=5, out='mag')
    current_image.show_image()

    current_image.contour_fitting()

    """
    image = current_image.image

    image_ga = gaussian(image, sigma=2)

    # seed = np.copy(image_ga)
    # seed[1:-1, 1:-1] = image_ga.min()

    image = (image - image.min()) / (image.max()-image.min())
    h = .9
    seed = image - h
    mask = image

    dilated = reconstruction(seed, mask, method='dilation')

    image = image - dilated

    global_thresh = threshold_otsu(image)
    binary_global = image > global_thresh

    block_size = 35
    adaptive_thresh = threshold_local(image, block_size, offset=10)
    binary_adaptive = image > adaptive_thresh

    binary_adaptive = np.invert(binary_adaptive)

    binary_adaptive = gaussian(binary_adaptive, sigma=1)

    fig, axes = plt.subplots(nrows=4, figsize=(7, 8))
    ax = axes.ravel()
    plt.gray()

    ax[0].imshow(image_ga)
    ax[0].set_title('Gaussian')

    ax[1].imshow(image)
    ax[1].set_title('Dilated')

    ax[2].imshow(binary_global)
    ax[2].set_title('Global thresholding')

    ax[3].imshow(binary_adaptive)
    ax[3].set_title('Adaptive thresholding')

    for a in ax:
        a.axis('off')

    plt.show()

    """
    """
    Contour detection with skimage
    im = (
        (current_image.image-current_image.image.min()) /
        (current_image.image.max() - current_image.image.min())
    )
    contours = measure.find_contours(im, 0.15)

    fig, ax = plt.subplots()
    ax.imshow(im, interpolation='nearest', cmap=plt.cm.gray)

    for n, contour in enumerate(contours):
        ax.plot(contour[:, 1], contour[:, 0], linewidth=2)

    ax.axis('image')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()
    """
    """
    current_image = ImageProcessing(path)
    current_image.equalize_hist_by_clahe(.2, 8)
    # current_image.gradient(type='scharr', size=5, out='mag')
    current_image.show_image()

    img = current_image.image
    img = np.array(img, dtype=np.uint8)
    entr_img = entropy(img, disk(30))
    global_thresh = threshold_otsu(entr_img)
    binary_global = entr_img > global_thresh

    fig, ax = plt.subplots(3, 1, figsize=(11, 7))
    ax[0].imshow(image, cmap='gray')

    img1 = ax[1].imshow(entr_img)
    ax[1].axis("off")

    ax[2].imshow(binary_global)

    fig.colorbar(img1, ax=ax[1])

    plt.show()
    """
