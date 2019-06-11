# -*- coding: utf-8 -*-
"""
Processing for solidification image analysis.

Authors
-------
Emeryk Ablonet : eablonet
    phD Student at IMFT[1]
    emeryk.ablonet@gmail.com
    [1]Institut de mecanique des fluides de Toulouse,
    2 allee du Professeur Camille Soula 31400 Toulouse

Copyright
---------
IMFT and Toulouse INsP

Dates
----
2018

Version
-------
8.1

Classes
-------
Stack

"""


# global packages
# ---------------
import sys
import os
import glob

import numpy as np

import progressbar

from matplotlib.backends.qt_compat import QtWidgets
from PySide2.QtWidgets import (QApplication, QFileDialog, QWidget)

from matplotlib import pyplot as plt

from scipy.optimize import curve_fit
from lmfit import Model

from skimage import morphology
from scipy.ndimage import morphology as morph
from skimage import measure
from scipy import ndimage

import pandas as ps

# local packages
# --------------
from . import ImageProcessing as ip
from . import read_online_data as rd

from ..geometrybuilder import LineBuilder as lb
from ..geometrybuilder import DraggrableRectangle as drag
from ..geometrybuilder import Geometry as geom

from ..gui import ClaheWindow as cw

# code
# ----


class Stack(object):
    """This class allow you to treat image by block.

    Parameters
    ----------
    data_directory : str
        directory containing images

    Methods
    -------
    lot of
    """

    def __init__(self):
        None

    def __str__(self):
        """Print the data directory on the terminal."""
        return self.data_directory

    def read_by_path(self, date, serie, regex='cam1_*', ext='.tif'):
        """Choose the path.

        If the experiement is not in the db, an error occur.

        Parameters
        ----------
        date: str
            select the date of the experiment
        serie : int
            select the serie of the experiment

        """
        self.data_directory = None
        df = rd.get_data()
        if date in list(df['date']):
            if serie in list(df[df['date'] == date]['serie']):
                self.data_directory = (
                    "/Volumes/EMERYK_HD/0_data/{:s}/n{:#}/".format(
                        date, int(serie)
                    )
                )
                self.datas = df[df['date'] == date][
                    df['serie'] == float(serie)
                ].copy()

        if self.data_directory is None:
            raise IndexError(
                'This experiment is not in the database, please add it before'
            )
        else:
            print('Path readed...')

        self.create_directory()
        self.update_lists(regex, ext)

        self.current_image_number = 1
        self.read_image()
        print('Current image selected : 1.')

    def select_folder(self):
        """
        UI select folder.

        todo
        ----

        """
        QApplication(sys.argv)
        folder = GetFolderName()
        serie = int(folder.name[-1])
        date = folder.name[-13:-3]

        self.read_by_path(date, serie)

    def create_directory(self):
        """Define output directories."""
        self.exporting_directory = self.data_directory + 'export_image/'

        if not os.path.isdir(self.exporting_directory):
            os.makedirs(self.exporting_directory + 'lastOp/')

    def update_lists(self, regex='cam1_*', ext='.tif'):
        """Update the lists of image, spatiox and Spatio_y."""
        # check what exists in data_folder #
        if os.path.isdir(self.data_directory + '_export'):
            # new procedure of reading # to work on pretreated images
            listdir = os.listdir(self.data_directory + '_export')
            print('\t' + str(0) + '. ' + 'Original reader image')
            c = 1
            for d in listdir:
                print('\t' + str(c) + '. ' + d)
                c += 1

            num = int(input('Select a folder to read : '))
            if num > 0:
                self.image_list = sorted(glob.glob((
                    self.data_directory + '_export/' + listdir[num-1] + '/' +
                    regex + ext
                )))
            else:
                self.image_list = sorted(glob.glob((
                    self.exporting_directory + 'lastOp/*.npy'
                )))
        else:  # all procedure
            # image list
            self.image_list = sorted(glob.glob((
                self.exporting_directory + 'lastOp/*.npy'
            )))

        if len(self.image_list) < 1:
            self.image_list = sorted(glob.glob((
                self.data_directory + regex + ext
            )))

        if num == 0:
            # image list
            self.image_list = sorted(glob.glob((
                self.exporting_directory + 'lastOp/*.npy'
            )))

            if len(self.image_list) < 1:
                self.image_list = sorted(glob.glob((
                    self.data_directory + regex + ext
                )))

        self.n_image_tot = len(self.image_list)
        self.range_images = range(1, self.n_image_tot+1)

    def print_info(self, type=None):
        """
        Show info about the stack.

        Print :
            1. directory which is loaded
            2. number of image loaded
            3. current image selected
            4. first image number of the stack
            5. last image number of the stack
        """
        print(
            '\nDirectory containing data :\n   {0}'.format(self.data_directory)
        )
        print(
            'Directory for exporting :\n   {0}'
            .format(self.exporting_directory)
        )
        print('---image info---')
        print('   Number of image : {0}'.format(self.n_image_tot))
        print(
            '   Current image selected : {0}'
            .format(self.current_image_number)
        )
        if type is not None:
            print('\n\n\n---List of images---')
            for i in range(len(self.image_list)):
                print('   #{0} : {1}'.format(i, self.image_list[i]))

        print('---end---\n\n')

    def read_image(self, ni=-1):
        """
        Read the image indicate by is number.

        Parameters
        ----------
        ni : int
            Image number to load. -1 to load current image image number.
            Default -1.

        """
        if ni == -1:
            ni = self.current_image_number

        elif ni in self.range_images:
            self.current_image_number = ni
        else:
            print(
                (
                    "ERROR reading image {}. It my be out of range\n" +
                    "Current image will be 1."
                ).format(
                    ni
                )
            )
            self.read_image(1)

        self.current_image = ip.ImageProcessing(
                self.image_list[ni-1]  # indexes starts at 0 in python
        )

    def treatment(self, treatment, *args, plot=True):
        """
        Apply a treatment to the current_image, and displays it.

        Only three treatments are available for now:
        ...crop : to crop the current_image (need 4 parameters:
                xmin, xmax, ymin, ymax)
        ...equalize : to equalize the image (global equalization, 2 parameters:
                limit = .2 by default, size= 8 by default)
        ...clahe : to equalize the image with clahe algorithm, local
                equalization, 2 parameters:
                    limit = .2 by default,
                    size= 8 by default
        """
        if treatment == 'crop':
            if len(args) == 4:
                self.current_image.crop(*args)
            else:
                raise StackEAError(
                    'Crop treatment requires 4 arguments'
                )
        elif treatment == 'clahe':
            if len(args) == 2 or len(args) == 0:
                self.current_image.equalize_hist_by_clahe(*args)
            else:
                raise StackEAError(
                    'Equalize treatment require 0 or 2 arguments'
                )
        elif treatment == 'equalize':
            if len(args) == 1 or len(args) == 0:
                self.current_image.equalize_hist(*args)
            else:
                raise StackEAError(
                    'Equalize treatment require 0 or 1 argument'
                )
        else:
            raise StackEAError(
                'The treatment you are trying to apply does not exist'
            )

        if plot:
            self.current_image.show_image()

    def remove_treatments(self):
        """Remove all lastOp files."""
        files = glob.glob((
            self.exporting_directory + 'lastOp/*.npy'
        ))
        for i in files:
            os.remove(i)
        f = open(self.data_directory + ".history", "a+")
        f.write(
            '\tRemove all files in {:s}\n'.format(
                self.exporting_directory + '/lastOp/'
            )
        )
        f.close()
        self.update_lists()
        self.read_image()

    def crop(self):
        """Crop with DraggableRectangle."""
        temp = self.current_image_number
        self.read_image(self.n_image_tot)

        fig, ax = plt.subplots(1)
        ax.imshow(self.current_image.image, cmap='gray')

        dr = drag.DraggableRectangle(ax)
        plt.show()

        if dr.valid:
            self.save_treatment('crop', *dr.get_rect())

        self.read_image(temp)

    def clahe(self):
        """Apply clahe thanks GUI."""
        qapp = QtWidgets.QApplication(sys.argv)
        app = cw.ClaheWindow(self.current_image)
        app.show()
        qapp.exec_()
        print(app.get_data())

        self.save_treatment('clahe', *app.get_data())

    def save_treatment(self, treatment, *args):
        """
        Treat all stakck image and save them in the results directory.

        Arguments are needed. There are the same than for treament, because
        save_treatment call recursivly treament.
        """
        temp = self.current_image_number

        # === initiate progressbar === #
        widgets = [treatment,
                   ' ', progressbar.Percentage(),
                   ' ', progressbar.Bar('=', '[', ']'),
                   ' ', progressbar.ETA(),
                   ' ', progressbar.FileTransferSpeed()]
        pbar = progressbar.ProgressBar(
            widgets=widgets, maxval=len(self.image_list)
        )
        pbar.start()
        # === loop over images === #
        for i in self.range_images:
            self.read_image(i)
            self.treatment(treatment, *args, plot=False)
            if i < 10:
                np.save(
                    self.exporting_directory + 'lastOp/' + 'i_000' + str(i),
                    self.current_image.image
                )
            elif i < 100:
                np.save(
                    self.exporting_directory + 'lastOp/' + 'i_00' + str(i),
                    self.current_image.image
                )
            elif i < 1000:
                np.save(
                    self.exporting_directory + 'lastOp/' + 'i_0' + str(i),
                    self.current_image.image
                )
            pbar.update(i)
        pbar.finish()

        # === update images infos === #
        self.update_lists()
        self.read_image(temp)

    def define_geometry(self):
        ref_pente = float(self.datas.iloc[0]['alpha'])
        temp = self.current_image_number
        self.read_image(0)
        im = self.current_image.rotate(angle=-ref_pente)
        self.read_image(temp)

        fig, ax = plt.subplots(1, figsize=[10, 7])
        ax.imshow(im, cmap='gray')

        dr = drag.DraggableRectangle(ax)
        plt.show()

        if dr.valid:
            rect = dr.get_rect()
            r0 = int((rect[1] - rect[0])/2)
            rc = int(r0 + rect[0])
            z0 = int(rect[3] - rect[2])

            print('\tz0 : ', z0)
            print('\tr0 : ', r0)
            print('\trc : ', rc)
            print('Please add this value in the table')

    def tracker(self, folder_name='front', pts=[]):
        """Get the front manually."""
        temp = self.current_image_number

        interpolation = np.zeros(
            (self.n_image_tot, self.current_image.size[1])
        )
        points = np.zeros(
            (self.n_image_tot, self.current_image.size[1])
        )
        r_space = range(0,  self.current_image.size[1])

        time_ref = int(self.datas.iloc[0]['t_ref_calc'])
        # time to start detection
        time_end = int(self.datas.iloc[0]['t_end_calc'])
        # time to end detection

        ref_pente = float(self.datas.iloc[0]['alpha'])
        # substrate inclination

        for n in range(time_ref, time_end):
            self.read_image(n)
            im = self.current_image.rotate(-ref_pente)

            # time gradient image ax3
            if n < self.n_image_tot+1 and n > 0:

                self.read_image(n+1)
                im_a = np.array(self.current_image.rotate(-ref_pente))
                self.read_image(n-1)
                im_b = np.array(self.current_image.rotate(-ref_pente))
                im_time = (2*im - im_b - im_a)/2

            elif n == 0:
                self.read_image(n+1)
                im_a = np.array(self.current_image.rotate(-ref_pente))
                im_time = im_a - im

            elif n == self.n_image_tot:
                self.read_image(n-1)
                im_b = np.array(self.current_image.rotate(-ref_pente))
                im_time = im - im_b

            im_grad = self.current_image.gradient(5, 'sobel', 'mag', 'valid')
            im_grady = self.current_image.gradient(5, 'sobel', 'y', 'valid')
            im_grad = self.current_image.rotate(-ref_pente, im_grad)

            fig = plt.figure(figsize=(20, 10))
            fig.canvas.set_window_title('Image {}/{}'.format(
                    n, self.n_image_tot+1
                )
            )
            ax1 = plt.subplot(2, 2, 1)
            ax1.imshow(im, cmap='gray')
            ax1.axis('off')
            ax1.set_title('Original image')

            ax2 = plt.subplot(2, 2, 2)
            ax2.imshow(im_grad, cmap='gray')
            ax2.axis('off')
            ax2.set_title('Gradient magnitude')

            ax3 = plt.subplot(2, 2, 3)
            ax3.imshow(im_time, cmap='gray')
            ax3.axis('off')
            ax3.set_title('Time gradient')

            ax4 = plt.subplot(2, 2, 4)
            ax4.imshow(im_grady, cmap='gray')
            ax4.axis('off')
            ax4.set_title('y-Gradient')

            if n > time_ref:
                # plot previous line
                ax1.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )
                ax2.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )
                ax3.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )
                ax4.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )

                # plot previous point
                ax1.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )
                ax2.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )
                ax3.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )
                ax4.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )

                ylim_low = np.min(linebuilder.ys[linebuilder.ys != 0]) - 40 # noqa: ignore=F821
                ylim_low = ylim_low if ylim_low > 0 else 0  # noqa: ignore=F821
                ylim_upp = np.max(linebuilder.ys[linebuilder.ys != 0]) + 10 # noqa: ignore=F821
                ylim_upp = ylim_upp if ylim_upp < im_grad.shape[0] else im_grad.shape[0]  # noqa: ignore=F821

                xlim_left = np.min(linebuilder.xs[linebuilder.xs != 0]) - 40 # noqa: ignore=F821
                xlim_left = xlim_left if xlim_left > 0 else 0  # noqa: ignore=F821
                xlim_right = np.max(linebuilder.xs[linebuilder.xs != 0]) + 40 # noqa: ignore=F821
                xlim_right = xlim_right if xlim_right < im_grad.shape[1] else im_grad.shape[1]  # noqa: ignore=F821

                ax2.set_xlim([xlim_left, xlim_right])
                ax2.set_ylim([ylim_upp, ylim_low])

                ax3.set_xlim([xlim_left, xlim_right])
                ax3.set_ylim([ylim_upp, ylim_low])

                ax4.set_xlim([xlim_left, xlim_right])
                ax4.set_ylim([ylim_upp, ylim_low])

            linebuilder = lb.SplineBuilder(ax1, ax2, ax3, ax4)
            plt.tight_layout()
            plt.show()

            points[n, linebuilder.xs] = linebuilder.ys
            interpolation[n, linebuilder.xs_interp] = linebuilder.ys_interp

        if not os.path.isdir(self.data_directory + folder_name):
            os.makedirs(self.data_directory + folder_name)

        if not os.path.isfile(self.data_directory + 'pente.npy'):
            np.save(
                self.data_directory + folder_name + '/interpolation',
                interpolation
            )
            np.save(
                self.data_directory + folder_name + '/points',
                points
            )
        else:
            answer = input(
                'Do you want to overwrite previous results ? [y/n] '
            )
            if answer.lower() in ['y', 'yes']:
                os.remove(self.data_directory + 'interpolation.npy')
                np.save(self.data_directory + 'interpolation', interpolation)
                os.remove(self.data_directory + 'points.npy')
                np.save(self.data_directory + 'points', points)

        self.read_image(temp)

        return points, interpolation

    def contour_tracker(self, folder_name='sa_contour', pts=[]):
        """Semi-automatic contour detection.

        Steps:
        -------
        1. Manual detect contour from tref to tref+4
        2. Manuel detect contour at (t_end-t_ref)/4, (t_end-t_ref)/2,
            3*(t_end-t_ref)/4, and t_end
        3. Get rc_left, rc_right, z0, and zf
        4. Create mask with previous data (M)
        5. Compute gradient for each image. Masked-it with M.
        6. Seperate left and right image.
        7. On each side get xc(yc) in px. (ndrl max du gradient dans le masque)
        8. Re-construct contour.

        """
        temp = self.current_image_number

        interpolation = np.zeros(
            (self.n_image_tot, self.current_image.size[1])
        )
        points = np.zeros(
            (self.n_image_tot, self.current_image.size[1])
        )
        Cr_left = np.zeros(
            (self.n_image_tot, self.current_image.size[0])
        )
        Cr_right = np.zeros(
            (self.n_image_tot, self.current_image.size[0])
        )
        Cz = np.zeros(
            (self.n_image_tot, self.current_image.size[1])
        )
        r_space = range(0,  self.current_image.size[1])

        time_ref = int(self.datas.iloc[0]['t_ref_calc'])
        # time to start detection
        time_end = int(self.datas.iloc[0]['t_end_calc'])
        time_nuc = int(self.datas.iloc[0]['t_nuc_calc'])
        # time to end detection

        ref_pente = float(self.datas.iloc[0]['alpha'])
        # substrate inclination

        im_list = [time_ref, time_nuc, time_end]
        # print(im_list)
        # im_list.append(int((time_end-time_ref)/4+time_ref))
        # im_list.append(int((time_end-time_ref)/2+time_ref))
        # im_list.append(int(3*(time_end-time_ref)/4+time_ref))
        # im_list.append(int(time_end))

        for n in im_list:
            self.read_image(n)
            im = self.current_image.rotate(-ref_pente)

            # time gradient image ax3
            if n < self.n_image_tot+1 and n > 0:

                self.read_image(n+1)
                im_a = np.array(self.current_image.rotate(-ref_pente))
                self.read_image(n-1)
                im_b = np.array(self.current_image.rotate(-ref_pente))
                im_time = (2*im - im_b - im_a)/2

            elif n == 0:
                self.read_image(n+1)
                im_a = np.array(self.current_image.rotate(-ref_pente))
                im_time = im_a - im

            elif n == self.n_image_tot:
                self.read_image(n-1)
                im_b = np.array(self.current_image.rotate(-ref_pente))
                im_time = im - im_b

            im_grad = self.current_image.gradient(5, 'sobel', 'mag', 'valid')
            im_grady = self.current_image.gradient(5, 'sobel', 'y', 'valid')
            im_grad = self.current_image.rotate(-ref_pente, im_grad)

            fig = plt.figure(figsize=(20, 10))
            fig.canvas.set_window_title('Image {}/{}'.format(
                    n, self.n_image_tot+1
                )
            )
            ax1 = plt.subplot(2, 2, 1)
            ax1.imshow(im, cmap='gray')
            ax1.axis('off')
            ax1.set_title('Original image')

            ax2 = plt.subplot(2, 2, 2)
            ax2.imshow(im_grad, cmap='gray')
            ax2.axis('off')
            ax2.set_title('Gradient magnitude')

            ax3 = plt.subplot(2, 2, 3)
            ax3.imshow(im_time, cmap='gray')
            ax3.axis('off')
            ax3.set_title('Time gradient')

            ax4 = plt.subplot(2, 2, 4)
            ax4.imshow(im_grady, cmap='gray')
            ax4.axis('off')
            ax4.set_title('y-Gradient')

            if n > time_ref:
                # plot previous line
                ax1.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )
                ax2.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )
                ax3.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )
                ax4.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )

                # plot previous point
                ax1.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )
                ax2.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )
                ax3.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )
                ax4.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )

                ylim_low = np.min(linebuilder.ys[linebuilder.ys != 0]) - 40 # noqa: ignore=F821
                ylim_low = ylim_low if ylim_low > 0 else 0  # noqa: ignore=F821
                ylim_upp = np.max(linebuilder.ys[linebuilder.ys != 0]) + 10 # noqa: ignore=F821
                ylim_upp = ylim_upp if ylim_upp < im_grad.shape[0] else im_grad.shape[0]  # noqa: ignore=F821

                xlim_left = np.min(linebuilder.xs[linebuilder.xs != 0]) - 40 # noqa: ignore=F821
                xlim_left = xlim_left if xlim_left > 0 else 0  # noqa: ignore=F821
                xlim_right = np.max(linebuilder.xs[linebuilder.xs != 0]) + 40 # noqa: ignore=F821
                xlim_right = xlim_right if xlim_right < im_grad.shape[1] else im_grad.shape[1]  # noqa: ignore=F821

                ax2.set_xlim([xlim_left, xlim_right])
                ax2.set_ylim([ylim_upp, ylim_low])

                ax3.set_xlim([xlim_left, xlim_right])
                ax3.set_ylim([ylim_upp, ylim_low])

                ax4.set_xlim([xlim_left, xlim_right])
                ax4.set_ylim([ylim_upp, ylim_low])

            linebuilder = lb.SplineBuilder(ax1, ax2, ax3, ax4)
            plt.tight_layout()
            plt.show()

            points[n, linebuilder.xs] = linebuilder.ys
            interpolation[n, linebuilder.xs_interp] = linebuilder.ys_interp

        # get geometry image time_ref
        fig = plt.figure(figsize=(20, 10))
        fig.canvas.set_window_title('Image {}/{}'.format(
                n, self.n_image_tot+1
            )
        )
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)

        self.read_image(time_ref)
        im = self.current_image.rotate(-ref_pente)
        ax1.imshow(im, cmap='gray')
        ax1.axis('off')

        self.read_image(time_end)
        im = self.current_image.rotate(-ref_pente)
        ax2.imshow(im, cmap='gray')
        ax2.axis('off')

        g = geom.Geometry(ax1, ax2, im.shape[0], im.shape[1])
        plt.tight_layout()
        plt.show()

        rc_l = g.rc_left_pos
        rc_r = g.rc_right_pos
        z0 = g.z0_pos
        zc = g.zc_pos
        print(rc_l, rc_r, z0, zc)

        rc_l = int(g.rc_left_pos)
        rc_r = int(g.rc_right_pos)
        rc = int(g.rc_pos)
        z0 = int(g.z0_pos)
        zc = int(g.zc_pos)
        zf = int(g.zf_pos)
        print(rc_l, rc_r, z0, zc)

        # create mask
        mask = im.copy()
        mask *= 0
        interpolation[interpolation == 0] = np.nan
        for r in range(rc_l, rc_r):
            if not np.isnan(np.nanmin(interpolation[:, r])):
                for z in range(
                    int(np.nanmin(interpolation[:, r]))-3,
                    int(np.nanmax(interpolation[:, r]))+3
                ):
                    mask[z, r] = 1

        plt.figure()
        plt.imshow(mask, cmap='gray')
        plt.show()

        for t in range(time_ref, time_end):
            self.read_image(t)
            im = self.current_image.rotate(-ref_pente)
            im_grad = self.current_image.gradient(5, 'sobel', 'mag', 'same')
            im_phi = self.current_image.gradient(5, 'sobel', 'arg', 'same')
            im_grad2 = self.current_image.gradient(5, 'sobel', 'mag', 'same', im_phi)

            # thres = np.median(im_grad) + .3*np.std(im_grad)

            # im_grad[im_grad < thres] = 0

            im_grad *= mask
            im_grad = np.ones_like(im_grad) - im_grad
            left = im_grad[:, :rc]
            right = im_grad[:, :rc:-1]
            if t == time_ref:
                plt.figure(figsize=[8, 4.5])
                plt.subplot(2,2,1)
                plt.imshow(left, cmap='gray')

                plt.subplot(2,2,2)
                plt.imshow(right, cmap='gray')

                plt.subplot(2,2,3)
                plt.imshow(im_phi*mask, cmap='gray')

                plt.subplot(2,2,4)
                plt.imshow(im_grad2*mask, cmap='gray')

            r_left = np.zeros(left.shape[0])
            r_right = np.zeros(left.shape[0])
            for z in range(0, left.shape[0]):
                c_l = np.argmin(left[z, :])
                if c_l > 0:
                    r_left[z] = c_l
                c_r = np.argmin(right[z, :])
                if c_r > 0:
                    r_right[z] = im_grad.shape[1] - c_r

            c_z = np.zeros(im_grad.shape[1])
            for r in range(rc_l, rc_r, 1):
                c = np.argmin(im_grad[:, r])
                if c > 0:
                    c_z[r] = c

            Cr_left[t, range(left.shape[0])] = r_left
            Cr_right[t, range(left.shape[0])] = r_right
            Cz[t, range(im_grad.shape[1])] = c_z

        self.read_image(time_ref+7)
        plt.figure(figsize=[8, 4.5])
        plt.imshow(self.current_image.image, cmap='gray')
        plt.plot(Cr_left[time_ref+7, :], range(0, left.shape[0]), '.c', ms=1)
        plt.plot(Cr_right[time_ref+7, :], range(0, right.shape[0]), '.m', ms=1)
        plt.plot(range(0, im_grad.shape[1]), Cz[time_ref+7, :], '.y', ms=1)
        plt.show()

        if not os.path.isdir(self.data_directory + folder_name):
            os.makedirs(self.data_directory + folder_name)

        if not os.path.isfile(self.data_directory + 'points.npy'):
            np.save(
                self.data_directory + folder_name + '/contour',
                Cr_left
            )
            np.save(
                self.data_directory + folder_name + '/interpolation',
                interpolation
            )
            np.save(
                self.data_directory + folder_name + '/points',
                points
            )
        else:
            answer = input(
                'Do you want to overwrite previous results ? [y/n] '
            )
            if answer.lower() in ['y', 'yes']:
                os.remove(
                    self.data_directory + folder_name + '/interpolation.npy'
                )
                np.save(
                    self.data_directory + folder_name + 'interpolation',
                    interpolation
                )
                os.remove(
                    self.data_directory + folder_name + 'points.npy'
                )
                np.save(self.data_directory + folder_name + 'points', points)
                os.remove(
                    self.data_directory + folder_name + 'contour.npy'
                )
                np.save(self.data_directory + folder_name + 'contour', Cr_left)

        self.read_image(temp)

        return Cr_left, points, interpolation

    def contour_tracker2(self, folder_name='sa_contour', pts=[]):
        """Semi-automatic contour detection.

        Steps:
        -------
        1. Get Geometry (rc_left, rc_right, rc, zc, z0, zf)
        2. Manual detect contour for time_nuc-1, time_nuc, time_end
            Preplace point of intersection rc_left, rc_right, z0 etc..
            Make a polyfit of 3rd order on left and right part.
        3.

        """
        temp = self.current_image_number

        interpolation = np.zeros(
            (self.n_image_tot, self.current_image.size[1])
        )
        points = np.zeros(
            (self.n_image_tot, self.current_image.size[1])
        )
        r_space = range(0,  self.current_image.size[1])

        time_ref = int(self.datas.iloc[0]['t_ref_calc'])
        # time to start detection
        time_end = int(self.datas.iloc[0]['t_end_calc'])
        time_nuc = int(self.datas.iloc[0]['t_nuc_calc'])
        # time to end detection

        ref_pente = float(self.datas.iloc[0]['alpha'])
        # substrate inclination

        """
        1. get geometry
        """
        fig = plt.figure(figsize=(20, 10))
        fig.canvas.set_window_title('time_ref and time_end images')
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)

        self.read_image(time_ref)
        im = self.current_image.rotate(-ref_pente)
        ax1.imshow(im, cmap='gray')
        ax1.axis('off')

        self.read_image(time_end)
        im = self.current_image.rotate(-ref_pente)
        ax2.imshow(im, cmap='gray')
        ax2.axis('off')

        g = geom.Geometry(ax1, ax2, im.shape[0], im.shape[1])
        plt.tight_layout()
        plt.show()

        rc_l = int(g.rc_left_pos)
        rc_r = int(g.rc_right_pos)
        rc = int(g.rc_pos)
        z0 = int(g.z0_pos)
        zc = int(g.zc_pos)
        zf = int(g.zf_pos)

        """
        2. Compute volume
        """
        im_list = [time_nuc-1, time_nuc, time_end]
        im_list = range(time_nuc-1, time_end)
        Cr_l = np.zeros(
            (self.n_image_tot, self.current_image.size[0])
        )
        Cr_r = np.zeros(
            (self.n_image_tot, self.current_image.size[0])
        )
        top_pt = np.zeros(self.n_image_tot)
        Volume = np.zeros(self.n_image_tot)
        for n in im_list:
            self.read_image(n)
            im = self.current_image.rotate(-ref_pente)
            im_grad = self.current_image.gradient(5, 'sobel', 'mag', 'valid')
            fig = plt.figure(figsize=(20, 10))
            fig.canvas.set_window_title('Image {}/{}'.format(
                    time_nuc, self.n_image_tot+1
                )
            )

            ax1 = plt.subplot(1, 2, 1)
            ax1.imshow(im, cmap='gray')
            ax1.axis('off')
            ax1.set_title('Original image')

            ax2 = plt.subplot(1, 2, 2)
            ax2.imshow(im_grad, cmap='gray')
            ax2.axis('off')
            ax2.set_title('Gradient magnitude')

            geom_data = [rc_l, rc_r, rc, zc, z0, zf]
            linebuilder = lb.ContourLineBuilder(ax1, ax2, geom_data)
            plt.tight_layout()
            plt.show()

            cr_l = np.array(linebuilder.ax1_left_line.get_xdata(), dtype=np.int)
            cz_l = np.array(linebuilder.ax1_left_line.get_ydata(), dtype=np.int)
            Cr_l[n, cz_l] = cr_l
            cr_r = np.array(linebuilder.ax1_right_line.get_xdata(), dtype=np.int)
            cz_r = np.array(linebuilder.ax1_right_line.get_ydata(), dtype=np.int)
            Cr_r[n, cz_r] = cr_r
            top_pt[n] = np.array(linebuilder.ax1_top_pt.get_ydata(), dtype=np.int)

            # volume calculation
            for i in range(self.current_image.size[0]):
                if Cr_l[n, i] != 0:
                    Volume[n] += (rc - Cr_l[n, i])**2
                if Cr_r[n, i] != 0:
                    Volume[n] += (Cr_r[n, i] - rc)**2

        dz_nuc = top_pt[time_nuc] - top_pt[time_nuc-1]
        dz_tot = top_pt[time_end] - top_pt[time_nuc-1]

        V0 = 0
        Vnuc = 0
        Vend = 0
        for i in range(self.current_image.size[0]):
            if Cr_l[time_nuc-1, i] != 0:
                V0 += (rc - Cr_l[time_nuc-1, i])**2
            if Cr_r[time_nuc-1, i] != 0:
                V0 += (Cr_r[time_nuc-1, i] - rc)**2

            if Cr_l[time_nuc, i] != 0:
                Vnuc += (rc - Cr_l[time_nuc, i])**2
            if Cr_r[time_nuc, i] != 0:
                Vnuc += (Cr_r[time_nuc, i] - rc)**2

            if Cr_l[time_end, i] != 0:
                Vend += (rc - Cr_l[time_end, i])**2
            if Cr_r[time_end, i] != 0:
                Vend += (Cr_r[time_end, i] - rc)**2

        V0 *= np.pi/2
        Vnuc *= np.pi/2
        Vend *= np.pi/2

        dV_nuc = Vnuc - V0
        dV_tot = Vend - V0

        print('dz_nuc : ', dz_nuc)
        print('dz_tot_geom', zf-z0)
        print('dz_tot : ', dz_tot)

        print('V0 :', V0)
        print('Vnuc :', Vnuc)
        print('Vend :', Vend)
        print('dV_nuc : ', dV_nuc)
        print('dV_tot : ', dV_tot)

        a_z = (998*top_pt[time_nuc-1]/top_pt[time_nuc] - 998) / (916-998)
        a_V = (998*V0/Vnuc - 998) / (916-998)

        print('alpha_z :', a_z)
        print('alpha_volume :', a_V)

        plt.figure(figsize=[8, 4.5])
        plt.plot(range(self.n_image_tot), Volume, '.b')
        # points[n, linebuilder.xs] = linebuilder.ys
        # interpolation[n, linebuilder.xs_interp] = linebuilder.ys_interp

        """
        for n in im_list:
            self.read_image(n)
            im = self.current_image.rotate(-ref_pente)

            # time gradient image ax3
            if n < self.n_image_tot+1 and n > 0:

                self.read_image(n+1)
                im_a = np.array(self.current_image.rotate(-ref_pente))
                self.read_image(n-1)
                im_b = np.array(self.current_image.rotate(-ref_pente))
                im_time = (2*im - im_b - im_a)/2

            elif n == 0:
                self.read_image(n+1)
                im_a = np.array(self.current_image.rotate(-ref_pente))
                im_time = im_a - im

            elif n == self.n_image_tot:
                self.read_image(n-1)
                im_b = np.array(self.current_image.rotate(-ref_pente))
                im_time = im - im_b

            im_grad = self.current_image.gradient(5, 'sobel', 'mag', 'valid')
            im_grady = self.current_image.gradient(5, 'sobel', 'y', 'valid')
            im_grad = self.current_image.rotate(-ref_pente, im_grad)

            fig = plt.figure(figsize=(20, 10))
            fig.canvas.set_window_title('Image {}/{}'.format(
                    n, self.n_image_tot+1
                )
            )
            ax1 = plt.subplot(2, 2, 1)
            ax1.imshow(im, cmap='gray')
            ax1.axis('off')
            ax1.set_title('Original image')

            ax2 = plt.subplot(2, 2, 2)
            ax2.imshow(im_grad, cmap='gray')
            ax2.axis('off')
            ax2.set_title('Gradient magnitude')

            ax3 = plt.subplot(2, 2, 3)
            ax3.imshow(im_time, cmap='gray')
            ax3.axis('off')
            ax3.set_title('Time gradient')

            ax4 = plt.subplot(2, 2, 4)
            ax4.imshow(im_grady, cmap='gray')
            ax4.axis('off')
            ax4.set_title('y-Gradient')

            if n > time_ref:
                # plot previous line
                ax1.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )
                ax2.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )
                ax3.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )
                ax4.plot(
                    r_space, interpolation[n-1, :],
                    '-r', alpha=.3
                )

                # plot previous point
                ax1.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )
                ax2.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )
                ax3.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )
                ax4.plot(
                    r_space, points[n-1, :],
                    'or', alpha=.3, ms=4,
                )

                ylim_low = np.min(linebuilder.ys[linebuilder.ys != 0]) - 40 # noqa: ignore=F821
                ylim_low = ylim_low if ylim_low > 0 else 0  # noqa: ignore=F821
                ylim_upp = np.max(linebuilder.ys[linebuilder.ys != 0]) + 10 # noqa: ignore=F821
                ylim_upp = ylim_upp if ylim_upp < im_grad.shape[0] else im_grad.shape[0]  # noqa: ignore=F821

                xlim_left = np.min(linebuilder.xs[linebuilder.xs != 0]) - 40 # noqa: ignore=F821
                xlim_left = xlim_left if xlim_left > 0 else 0  # noqa: ignore=F821
                xlim_right = np.max(linebuilder.xs[linebuilder.xs != 0]) + 40 # noqa: ignore=F821
                xlim_right = xlim_right if xlim_right < im_grad.shape[1] else im_grad.shape[1]  # noqa: ignore=F821

                ax2.set_xlim([xlim_left, xlim_right])
                ax2.set_ylim([ylim_upp, ylim_low])

                ax3.set_xlim([xlim_left, xlim_right])
                ax3.set_ylim([ylim_upp, ylim_low])

                ax4.set_xlim([xlim_left, xlim_right])
                ax4.set_ylim([ylim_upp, ylim_low])

            linebuilder = lb.SplineBuilder(ax1, ax2, ax3, ax4)
            plt.tight_layout()
            plt.show()

            points[n, linebuilder.xs] = linebuilder.ys
            interpolation[n, linebuilder.xs_interp] = linebuilder.ys_interp
        """

        """ old contour detection

        # get geometry image time_ref
        fig = plt.figure(figsize=(20, 10))
        fig.canvas.set_window_title('Image {}/{}'.format(
                n, self.n_image_tot+1
            )
        )
        ax1 = plt.subplot(1, 2, 1)
        ax2 = plt.subplot(1, 2, 2)

        self.read_image(time_ref)
        im = self.current_image.rotate(-ref_pente)
        ax1.imshow(im, cmap='gray')
        ax1.axis('off')

        self.read_image(time_end)
        im = self.current_image.rotate(-ref_pente)
        ax2.imshow(im, cmap='gray')
        ax2.axis('off')

        g = geom.Geometry(ax1, ax2, im.shape[0], im.shape[1])
        plt.tight_layout()
        plt.show()

        rc_l = int(g.rc_left_pos)
        rc_r = int(g.rc_right_pos)
        rc = int(g.rc_pos)
        z0 = int(g.z0_pos)
        zc = int(g.zc_pos)
        zf = int(g.zf_pos)
        print(rc_l, rc_r, z0, zc)

        # create mask
        mask = im.copy()
        mask *= 0
        interpolation[interpolation == 0] = np.nan
        for r in range(rc_l, rc_r):
            if not np.isnan(np.nanmin(interpolation[:, r])):
                for z in range(
                    int(np.nanmin(interpolation[:, r]))-3,
                    int(np.nanmax(interpolation[:, r]))+3
                ):
                    mask[z, r] = 1

        plt.figure()
        plt.imshow(mask, cmap='gray')
        plt.show()

        for t in range(time_ref, time_end):
            self.read_image(t)
            im = self.current_image.rotate(-ref_pente)
            im_grad = self.current_image.gradient(5, 'sobel', 'mag', 'same')

            # thres = np.median(im_grad) + .3*np.std(im_grad)

            # im_grad[im_grad < thres] = 0

            im_grad *= mask

            left = im_grad[:, :rc]
            right = im_grad[:, rc:]
            if t == time_ref:
                plt.figure(figsize=[8, 4.5])
                plt.subplot(1,2,1)
                plt.imshow(left)

                plt.subplot(1,2,2)
                plt.imshow(right)

            r_left = np.zeros(left.shape[0])
            r_right = np.zeros(left.shape[0])
            for z in range(0, left.shape[0]):
                c_l = np.argmax(left[z, :])
                if c_l > 0:
                    r_left[z] = c_l
                c_r = np.argmax(right[z, :])
                if c_r > 0:
                    r_right[z] = rc + c_r

            Cr_left[t, range(left.shape[0])] = r_left
            Cr_right[t, range(left.shape[0])] = r_right

        self.read_image(time_ref+7)
        plt.figure(figsize=[8, 4.5])
        plt.imshow(self.current_image.image, cmap='gray')
        plt.plot(Cr_left[time_ref+7, :], range(0, left.shape[0]), '.c', ms=1)
        plt.plot(Cr_right[time_ref+7, :], range(0, right.shape[0]), '.m', ms=1)
        plt.show()
        """

        if not os.path.isdir(self.data_directory + folder_name):
            os.makedirs(self.data_directory + folder_name)

        if not os.path.isfile(self.data_directory + 'volume.npy'):
            np.save(
                self.data_directory + folder_name + '/volume',
                Volume
            )
            np.save(
                self.data_directory + folder_name + '/Cr_left',
                Cr_l
            )
            np.save(
                self.data_directory + folder_name + '/Cr_right',
                Cr_r
            )
        else:
            answer = input(
                'Do you want to overwrite previous results ? [y/n] '
            )
            if answer.lower() in ['y', 'yes']:
                os.remove(
                    self.data_directory + folder_name + '/interpolation.npy'
                )
                np.save(
                    self.data_directory + folder_name + 'interpolation',
                    interpolation
                )
                os.remove(
                    self.data_directory + folder_name + 'points.npy'
                )
                np.save(self.data_directory + folder_name + 'points', points)
                os.remove(
                    self.data_directory + folder_name + 'contour.npy'
                )
                np.save(self.data_directory + folder_name + 'contour', Cr_r)

        self.read_image(temp)

        # return Cr_left, points, interpolation

    def view_all_profil(self, zf, n_space):
        """Plot front dynamic and propagation."""

        px_mm = float(self.datas.iloc[0]['px_mm'])
        t_ref = int(self.datas.iloc[0]['t_ref_calc'])
        rc = int(self.datas.iloc[0]['rc'])

        zf[zf == 0] = np.nan
        zf = self.current_image.size[0] - zf

        zf -= np.nanmean(zf[t_ref, :])

        fig, ax0 = plt.subplots(figsize=[8, 6])
        fig, ax1 = plt.subplots(figsize=[8, 6])
        x = np.arange(self.current_image.size[1])
        x = x[::n_space] - rc
        # colormap = plt.get_cmap('spring')
        # ax1.set_prop_cycle(
        #     [
        #         colormap(i) for i in np.linspace(0, 0.9, len(zf[:, 0]))
        #     ]
        # )

        for time in range(len(zf[:, 0])):
            ax1.plot(
                x/px_mm*1e3, zf[time, ::n_space]/px_mm*1e3,
                '--', linewidth=2,
            )
            ax0.plot(
                x, zf[time, ::n_space],
                '--', linewidth=2
            )
        ax0.set_xlabel('r (px)')
        ax0.set_ylabel('z (px)')
        ax0.grid(True)

        ax1.set_xlabel('r (mm)', fontsize=18)
        ax1.set_ylabel('z (mm)', fontsize=18)
        ax1.grid(True)
        plt.show()

    def read_manual_front(self):
        """Plot front dynamic and propagation.

        Returns
        slope : ndarray
            Mean inclination of the the solidificaito front at each time
        zf : ndarray
            Position of the solidification front at each time and for each x
            position

        """
        slope = np.load(self.data_directory + 'pente.npy')
        zf = np.load(self.data_directory + 'y0.npy')
        zf[zf == 0] = np.nan

        return slope, zf

    def read_data(self, folder_name='front'):
        """Plot front dynamic and propagation.

        Returns
        zf : ndarray
            Position of the solidification front at each time and for each x
            position
        """
        zf = np.load(self.data_directory + folder_name + '/interpolation.npy')
        # pts = np.load(self.data_directory + folder_name + '/points.npy')
        zf[zf == 0] = np.nan

        return zf#, pts

    def read_sa_contour(self, folder_name='sa_contour'):
        """Plot front dynamic and propagation.

        Returns
        zf : ndarray
            Position of the solidification front at each time and for each x
            position
        """
        vol = np.load(self.data_directory + folder_name + '/volume.npy')
        Cl = np.load(self.data_directory + folder_name + '/Cr_left.npy')
        Cr = np.load(self.data_directory + folder_name + '/Cr_right.npy')
        # pts = np.load(self.data_directory + folder_name + '/points.npy')
        # zf[zf == 0] = np.nan

        return vol, Cl, Cr  #, pts

    def threshold(self):
        self.current_image.plot_hist()
        plt.show()
        a = int(input('Select low threshold value : \n'))
        b = int(input('Select upper threshold value : \n'))

        # === loop over images === #
        for i in range(self.n_image_tot):
            self.read_image(i)
            im = self.current_image.image
            im[im < a] = a
            im[im > b] = b

            if i < 9:
                np.save(
                    self.exporting_directory + 'lastOp/' + 'i_000' + str(i+1),
                    im
                )
            elif i < 99:
                np.save(
                    self.exporting_directory + 'lastOp/' + 'i_00' + str(i+1),
                    im
                )
            elif i < 999:
                np.save(
                    self.exporting_directory + 'lastOp/' + 'i_0' + str(i+1),
                    im
                )

    def hist_evol(self):

        fig, [ax1, ax2] = plt.subplots(nrows=2, ncols=1, figsize=[11, 7])
        for i in range(self.n_image_tot):
            self.read_image(i)
            H, cs = self.current_image.plot_hist()
            print(np.mean(H))
            ax1.plot(i, np.mean(H), '.r')
            ax2.plot(i, np.std(H), '.r')
        plt.show()

    def make_montage(
        self, images, add_front=False, add_timer=False, rot=False
    ):
        """Create a montage.

        Parameters
        ----------
        image : list of int
            List of all image to plot.

        """
        rows = (len(images)-1) // 4 + 1
        if rows > 1:
            cols = 4
        else:
            cols = len(images) % 4

        fig = plt.figure(
            figsize=(2*cols, 2*rows)
        )
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None)

        if add_front:
            y0 = np.load(self.data_directory + 'y0.npy')
            y0[y0 == 0] = np.nan
            x_rot = np.arange(self.current_image.size[1])
            ref_pente = self.get_data(['alpha'])[0]

        if add_timer:
            time_ref, fps = self.get_data(
                ['time_ref', 'fps']
            )

        temp = self.current_image_number
        for row in range(rows):
            for col in range(cols):
                k = row*cols + col
                if k < len(images):
                    self.read_image(images[k])
                    if rot:
                        im = self.current_image.rotate(-ref_pente)
                    else:
                        im = self.current_image.image

                    ax = plt.subplot(rows, cols, k+1)
                    ax.imshow(im, cmap='gray')
                    if add_front:
                        if rot:
                            x = x_rot
                            y = y0[images[k], :]
                        else:
                            x, y = self._rot_inv(
                                x_rot,
                                y0[images[k], :],
                                -ref_pente
                            )
                        ax.plot(
                            x, y,
                            '.g', markersize=1,
                        )
                    if add_timer:
                        ax.text(
                            5,
                            5,
                            '{:.1f}'.format((images[k] - time_ref-4)/fps),
                            bbox={'facecolor': 'white', 'pad': 5}
                        )
                    ax.axis('off')

        fig.tight_layout()
        # gridspec_kw={'wspace': .05, 'hspace': .05}, squeeze=True
        plt.show()
        self.read_image(temp)

    def view_propagation_front(self, zf):
        """Display front inclination over timeself."""

        """
        time_ref, fps, ref_pente = self.get_data(
            ['time_ref', 'fps', 'alpha']
        )
        """

        time_ref = int(self.datas['t_ref_calc'].values)
        fps = float(self.datas.iloc[0]['fps'].replace(',', '.'))
        ref_pente = float(self.datas.iloc[0]['alpha'])

        zf = zf[time_ref:, :]

        def fit_fun(x, a, b):
            zf = a*x + b

            return zf

        p = np.zeros(len(zf[:, 0]))
        perr = np.zeros(len(zf[:, 0]))

        for i in range(len(zf[:, 0])):
            idx = np.isfinite(zf[i, :])
            popt, pcov = curve_fit(
                fit_fun,
                np.arange(0, len(zf[i, idx])),
                zf[i, idx]
            )
            p[i] = popt[0]
            perr[i] = np.sqrt(np.diag(pcov))[0]

        f, ax = plt.subplots(figsize=[11, 7])
        plt.errorbar(
            np.arange(len(p))/fps,
            np.arctan(p)*180/np.pi + ref_pente,
            yerr=np.arctan(perr)*180/np.pi,
            marker='s', markersize=4, markeredgecolor='k',
            markerfacecolor='None',
            linestyle='None',
            ecolor='k', elinewidth=.5
        )

        plt.plot(
            [0, len(p)/fps],
            [
                np.nanmean(np.arctan(p))+ref_pente,
                np.nanmean(np.arctan(p))+ref_pente
            ],
            '--k'
        )
        plt.plot(
            [0, len(p[time_ref:])/fps],
            [ref_pente, ref_pente],
            '--k', linewidth=2
        )
        plt.xlabel('Time (s)', fontsize=14)
        plt.ylabel('Inclination (°)', fontsize=14)
        plt.grid(True)

        ni = int(self.n_image_tot/4)
        temp_ni = self.current_image_number
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        xlim = list(xlim)
        ylim = list(ylim)
        xlim[0] += .05*xlim[0] if xlim[0] > 0 else -xlim[0]
        xlim[1] -= .05*xlim[1]
        ylim[1] -= .05*ylim[1]

        """
        self.read_image(ni)
        ax_ins = inset_axes(ax, width="30%", height="30%", loc=2)
        ax_ins.imshow(self.current_image.image, cmap='gray')
        plt.axis('off')

        ax.annotate(
            '', xy=(ni/fps, np.arctan(p[ni])*180/np.pi), xytext=(xlim[0], ylim[1]),
            arrowprops=dict(facecolor='black', arrowstyle="->"),
        )

        self.read_image(2*ni)
        ax_ins = inset_axes(ax, width="30%", height="30%", loc=9)
        ax_ins.imshow(self.current_image.image, cmap='gray')
        plt.axis('off')

        ax.annotate(
            '', xy=(2*ni/fps, np.arctan(p[2*ni])*180/np.pi), xytext=((xlim[1]-xlim[0])/2, ylim[1]),
            arrowprops=dict(facecolor='black', arrowstyle="->"),
        )

        self.read_image(3*ni)
        ax_ins = inset_axes(ax, width="30%", height="30%", loc=1)
        ax_ins.imshow(self.current_image.image, cmap='gray')
        plt.axis('off')

        ax.annotate(
            '', xy=(3*ni/fps, np.arctan(p[3*ni])*180/np.pi), xytext=(xlim[1], ylim[1]),
            arrowprops=dict(facecolor='black', arrowstyle="->"),
        )
        """

        plt.show()

        self.read_image(temp_ni)

    def get_dynamic_front(self, zf, usi=False):
        """Return correct dynamic front.

        Parameters
        ----------
        zf : ndarray
            Contain "brut" front information
        usi : bool
            False by default. If true, convert pixel length into meter,
            and time in second

        """
        time_ref = int(self.datas['t_ref_calc'].values)
        px_mm = float(self.datas.iloc[0]['px_mm'])
        fps = float(self.datas.iloc[0]['fps'])
        x_space = np.arange(self.current_image.size[1])
        t_space = np.arange(len(zf[time_ref:, 0]))

        zf = zf[time_ref:, :]
        zf_mean = np.zeros(len(zf[:, 0]))
        zf_std = np.zeros(len(zf[:, 0]))

        for i in range(1, len(zf[:, 0])):
            for j in x_space:
                if zf[i, j] != 0:
                    zf[i, j] = (zf[0, j] - zf[i, j])
                else:
                    zf[i, j] = np.nan
            if all(np.isnan(zf[i, :])):  # if there is only nan
                zf_mean[i] = np.inf
                zf_std[i] = np.nan
            else:
                zf_mean[i] = np.nanmean(zf[i, :])
                zf_std[i] = np.nanstd(zf[i, :])
        zf[0, :] = 0
        zf_mean[0] = 0

        # zf_mean[zf_mean == np.inf] = np.nan

        if usi:
            t_space = t_space/fps
            zf /= px_mm
            zf_mean /= px_mm
            zf_std /= px_mm

        return t_space, zf, zf_mean, zf_std

    def get_dynamic_front_by_loc(self, zf, loc, usi=False):
        """Return correct dynamic front for a specific location.

        Parameters
        ----------
        zf : ndarray
            Contain "brut" front information
        loc : int
            Position in e_r to return (in pixel !)
        usi : bool
            False by default. If true, convert pixel length into meter,
            and time in second

        """
        time_ref = int(self.datas['t_ref_calc'].values)
        px_mm = float(self.datas.iloc[0]['px_mm'])

        zf = zf[time_ref:, :]
        zf_loc = np.zeros(len(zf[:, 0]))

        for i in range(1, len(zf[:, 0])):
            if zf[i, loc] != 0:
                zf_loc[i] = (zf[0, loc] - zf[i, loc])
            else:
                zf_loc[i] = np.nan

        zf_loc[0] = 0

        zf_loc[zf_loc == np.inf] = np.nan

        if usi:
            zf_loc /= px_mm

        return zf_loc

    def get_dynamic_front_by_win(self, zf, center, width, usi=False):
        """Return correct dynamic front in a window.

        Parameters
        ----------
        zf : ndarray
            Contain "brut" front information
        center : int
            Position of center in e_r to return.
        width : int
            +- pixel to add on left and right of center
        usi : bool
            False by default. If true, convert pixel length into meter,
            and time in second

        """
        time_ref = int(self.datas['t_ref_calc'].values)
        px_mm = float(self.datas.iloc[0]['px_mm'])

        zf = zf[time_ref:, :]  # we don't care before time_ref

        r_space = range(center-width, center+width+1)
        zf_win = np.zeros([len(zf[:, 0]), len(r_space)])
        zf_win_mean = np.zeros(len(zf[:, 0]))
        zf_win_std = np.zeros(len(zf[:, 0]))

        for i in range(1, len(zf[:, 0])):
            for id, val in enumerate(r_space):
                if zf[i, val] != 0:
                    zf_win[i, id] = (zf[0, val] - zf[i, val])
                else:
                    zf_win[i, id] = np.nan

            zf_win_mean[i] = np.nanmean(zf_win[i, :])
            zf_win_std[i] = np.nanstd(zf_win[i, :])

        zf_win_mean[0] = 0
        zf_win_mean[zf_win_mean == np.inf] = np.nan

        if usi:
            zf_win_mean /= px_mm
            zf_win_std /= px_mm

        return zf_win_mean, zf_win_std

    def get_fit_front(self, t, zf, rho_s, Lf, k_s, fit_win=[5, 10]):
        delta = []
        t0 = []
        dT = []

        # --- fit with lmfit.model --- #
        def fit_fun_dir(t, delta, t0):
            zf = delta**2*(t - t0)
            return zf

        gmodel = Model(fit_fun_dir)
        gmodel.eval(
            t=t[fit_win[0]:-fit_win[1]],
            delta=1e-4, t0=.1,
        )
        result = gmodel.fit(
            (zf[fit_win[0]:-fit_win[1]])**2,
            t=t[fit_win[0]:-fit_win[1]],
            delta=1e-4, t0=.1,
        )
        delta.append(result.best_values['delta'])
        t0.append(result.best_values['t0'])
        dT.append(delta[0]**2*rho_s*Lf/(2*k_s))

        # --- fit with polyfit --- #
        p = np.polyfit(
            t[fit_win[0]:-fit_win[1]],
            (zf[fit_win[0]:-fit_win[1]]*1e3)**2,
            1
        )
        delta.append(np.sqrt(p[0]*1e-6))
        t0.append(-p[1]*1e-6/p[0]**2)
        dT.append(p[0]**2*rho_s*Lf/(2*k_s))

        return delta, t0, dT

    def get_fits(self, t, zf, plot=False):
        """
            Return to fit coefficient for t0 t0-1
        """
        t0, time_ref, fps = self.get_data(['t0', 'time_ref', 'fps'])
        t0 -= 1
        end_time = int(round(len(t)*.15))

        print(zf)
        print(t0-time_ref)
        print(end_time)

        def fit_fun_dir(t, Delta):
            zf_square = Delta*(t)
            return zf_square

        gmodel = Model(fit_fun_dir)
        gmodel.eval(
            t=t[t0-time_ref:-end_time]-(t0-time_ref)/fps,
            Delta=1e-4,
        )
        result = gmodel.fit(
            (zf[t0-time_ref:-end_time])**2,
            t=t[t0-time_ref:-end_time]-(t0-time_ref)/fps,
            Delta=1e-4,
        )
        delta_fit = result.best_values['Delta']
        delta_fit = np.sqrt(delta_fit)
        print(delta_fit)

        gmodel = Model(fit_fun_dir)
        gmodel.eval(
            t=t[t0-time_ref:-end_time]-(t0-1-time_ref)/fps,
            Delta=1e-4,
        )
        result = gmodel.fit(
            (zf[t0-time_ref:-end_time])**2,
            t=t[t0-time_ref:-end_time]-(t0-1-time_ref)/fps,
            Delta=1e-4,
        )
        delta_fit_down = result.best_values['Delta']
        delta_fit_down = np.sqrt(delta_fit_down)
        print(delta_fit_down)

        if plot:
            fig, ax0 = plt.subplots(figsize=(20, 10))
            ax0.plot(t-(t0-time_ref)/fps, zf**2, 'ok', markersize=8, markerfacecolor='None')
            ax0.plot(
                t[t0-time_ref:-end_time]-(t0-time_ref)/fps,
                zf[t0-time_ref:-end_time]**2,
                '*b',
            )
            ax0.plot(
                t, delta_fit**2*(t-(t0-time_ref)/fps), '--b'
            )
            ax0.plot(
                t[t0-time_ref:-end_time]-(t0-1-time_ref)/fps,
                (zf[t0-time_ref:-end_time])**2,
                '*r',
            )
            ax0.plot(
                t, delta_fit_down**2*(t-(t0-1-time_ref)/fps), '--r'
            )

            fig, ax = plt.subplots(figsize=(20, 10))
            ax.plot(t, zf, '*k')
            ax.plot(t, delta_fit*np.sqrt(t-(t0-time_ref)/fps), '--b')
            ax.plot(t, delta_fit_down*np.sqrt(t-(t0-1-time_ref)/fps), '--r')

            plt.show()

    def get_limit_stefan(self, t, zf, plot=False):
        """
            Return curves between t0 & t0+1 of Stefan profile.
        """
        t0, time_ref, fps, Tc, Tcs = self.get_data(['t0', 'time_ref', 'fps', 'Tc', 'Tcs'])
        t0 -= 1
        end_time = int(round(len(t)*.15))

        print(zf)
        print(t0-time_ref)
        print(end_time)

        # theoritical plot #
        # ---------------- #
        liquid = ste2.Material()
        liquid.set_data([998, 4200, 0.5])
        ice = ste2.Material()
        ice.set_data([910, 2050, 2.22])
        liquid_ice = ste2.PhaseTransitonMaterial(liquid, ice)
        liquid_ice.set_melting_temperature(0)
        liquid_ice.set_heat_latent(335000)
        th = ste2.Stefan(liquid_ice)
        zf_Tc = th.monophasic_steady(Tc, z0=2.2e-3)
        if Tcs is not None:
            zf_Tcs = th.monophasic_steady(Tcs, z0=2.2e-3)

        plt.plot(
            np.arange(len(zf))*th.dt + t0, zf,
            '--k',
            label=(
                r'Theoretical front $\delta$ = {:.1e}'.format(
                    np.sqrt(2*-Tc*k_s/(rho_s*Lf))
                )
            )
        )

    def view_dynamic_front(self, zf, fit_win=[0, 10]):
        """Display front dynamic.

        Parameters
        ----------
        zf : ndarray
            Contain front information
        fit_win : list of 2 elements
            Windowing the elements for fit. First element is the number of
            point to exlude at begin, and second is the number to exclude at
            the end.
        plot : bool
            True by default. Enable to plot results.

        """
        time_ref = int(self.datas['t_nuc_calc'].values)
        px_mm = int(self.datas['px_mm'].values)
        fps = float(self.datas.iloc[0]['fps'])
        try:
            Tc = float(self.datas.iloc[0]['Tc_nuc'])
        except AttributeError:
            Tc = float(self.datas.iloc[0]['Tc_set'])

        ref_pente = float(self.datas.iloc[0]['alpha'])

        z0 = self.datas['z0'].values
        r0 = self.datas['r0'].values
        rc = self.datas['rc'].values

        try:
            z0 = int(z0)
            r0 = int(r0)
            rc = int(rc)
        except ValueError:

            temp = self.current_image_number
            self.read_image(0)
            im = self.current_image.rotate(angle=-ref_pente)
            self.read_image(temp)

            fig, ax = plt.subplots(1, figsize=[10, 7])
            ax.imshow(im, cmap='gray')

            dr = drag.DraggableRectangle(ax)
            plt.show()

            if dr.valid:
                rect = dr.get_rect()
                r0 = int((rect[1] - rect[0])/2)
                rc = int(r0 + rect[0])
                z0 = int(rect[3] - rect[2])

                print('\tz0 : ', z0)
                print('\tr0 : ', r0)
                print('\trc : ', rc)

            self.read_image(temp)

        rho_s = 916.8
        Lf = 334000.
        k_s = 2.22

        t_space, zf, zf_mean, zf_std = self.get_dynamic_front(zf)

        t0idx = np.argmin(abs(t_space/fps-fit_win[0]))

        # plot for different rows #
        # ----------------------- #
        cols = np.array(
            np.linspace(rc-r0, rc+r0, 6, endpoint=False),
            dtype=np.int
        )
        cols = np.delete(cols, 0)
        plt.figure(figsize=[8, 4.5])
        for col in cols:
            plt.plot(
                np.array(t_space)/np.max(t_space),
                zf[:, col]/z0,
                's', markerfacecolor='None',
                label=r'$r/r_0 =$ {0:.2f}'.format((col-rc)/r0)
            )
        plt.xlabel(r'$t/t_f$')
        plt.ylabel(r'$z_f(t)/z_0$')
        plt.grid(True)
        plt.legend(loc=0, fancybox=True)

        plt.figure(figsize=[8, 4.5])
        self.read_image(time_ref+30)
        im = self.current_image.image
        plt.imshow(self.current_image.rotate(-ref_pente, im), cmap='gray')
        for col in cols:
            plt.plot(
                [col, col],
                [0, self.current_image.size[0]],
                '-.', markersize=2, alpha=.3,
            )

        """
        # model and fit comparison #
        # ------------------------ #
        # --- fit direct --- #
        def fit_fun_dir(t, delta, t0):
            zf = delta**2*(t - t0)
            return zf

        gmodel = Model(fit_fun_dir)
        print(zf_mean[fit_win[0]:-fit_win[1]]/px_mm)
        print(t_space[fit_win[0]:-fit_win[1]]/fps)
        gmodel.eval(
            t=t_space[fit_win[0]:-fit_win[1]]/fps,
            delta=1e-4, t0=.1,
        )
        result = gmodel.fit(
            (zf_mean[fit_win[0]:-fit_win[1]]/px_mm)**2,
            t=t_space[fit_win[0]:-fit_win[1]]/fps,
            delta=1e-4, t0=.1,
        )
        print(result.fit_report())
        print(result.best_values)
        delta = result.best_values['delta']
        t0 = result.best_values['t0']
        dT_dir = delta**2*rho_s*Lf/(2*k_s)

        # --- fit **2 -- #
        def fit_fun(t, delta, t0):
            zf_square = delta**2*(t - t0)
            return zf_square

        popt, pcov = curve_fit(
            fit_fun,
            t_space[fit_win[0]:-fit_win[1]]/fps,
            (zf_mean[fit_win[0]:-fit_win[1]]/px_mm*1e3)**2,
            bounds=(
                [0, 0],
                [np.inf, 6/fps]
            ),
            p0=[1e-4, t_space[fit_win[0]]/fps - 1/(2*fps)],
        )
        popt[0] = popt[0]*1e-3

        dT_square = popt[0]**2*rho_s*Lf/(2*k_s)
        perr = np.sqrt(np.diag(pcov))
        print(
            [0, t_space[fit_win[0]]/fps-fps],
            [np.inf, t_space[fit_win[0]]/fps]
        )
        print('Error :{}'.format(perr))

        p = np.polyfit(
            t_space[fit_win[0]:-fit_win[1]]/fps,
            (zf_mean[fit_win[0]:-fit_win[1]]/px_mm*1e3)**2,
            1
        )
        p[0] = np.sqrt(p[0]*1e-6)
        p[1] = -p[1]*1e-6/p[0]**2
        dT_poly = p[0]**2*rho_s*Lf/(2*k_s)
        """

        # -- plot fig --- #
        plt.figure(figsize=(8, 4.5))
        plt.errorbar(
            t_space/fps - 0/fps,
            zf_mean/px_mm,
            yerr=(zf_std+5)/px_mm,
            linestyle="None",
            marker='s',
            markersize=4, markerfacecolor='None'
        )
        plt.plot(
            t_space[t0idx:-fit_win[1]]/fps - 0/fps,
            zf_mean[t0idx:-fit_win[1]]/px_mm,
            '.r', markersize=4,
            label='Considered point for fitting'
        )

        z = np.polyfit(t_space[t0idx:-fit_win[1]]/fps, zf_mean[t0idx:-fit_win[1]]/px_mm, 1)
        p = np.poly1d(z)
        plt.plot(t_space/fps, p(t_space/fps), '--k', label='qw = '+str(z[0]*rho_s*Lf))

        """
        # --- fit dir --- #
        plt.plot(
            t_space[fit_win[0]:-fit_win[1]]/fps,
            delta*np.sqrt(t_space[fit_win[0]:-fit_win[1]]/fps-t0),
            'r-',
            label=r'lmodel $z_f^2 = \delta^2(t-t_0)$, $\delta$ = {:.1e}, $t_0$ = {:.1e}, $\Delta T$ = {:.1e}'.format(
                delta, t0, dT_dir
            )
        )
        """

        """
        # --- fit **2 ---#
        plt.plot(
            t_space/fps,
            popt[0]*np.sqrt(t_space/fps - popt[1]),
            '-.r',
            label=r'curve_fit $z_f^2 = \delta^2(t-t_0)$, $\delta$ = {:.1e}, $t_0$ = {:.1e}, $\Delta T$ = {:.1e}'.format(
                popt[0], popt[1], dT_square
            )
        )
        """

        """
        plt.plot(
            t_space/fps,
            p[0]*np.sqrt(t_space/fps - p[1]),
            '-.m',
            label=r'polyfit $z_f^2 = \delta^2*t - z0$, $\delta$ = {:.1e}, $t_0$ = {:.1e}, $\Delta T$ = {:.1e}'.format(
                p[0], p[1], dT_poly
            )
        )
        """

        """
        # theoritical plot #
        # ---------------- #
        liquid = ste2.Material()
        liquid.set_data([998, 4200, 0.5])
        ice = ste2.Material()
        ice.set_data([910, 2050, 2.22])
        liquid_ice = ste2.PhaseTransitonMaterial(liquid, ice)
        liquid_ice.set_melting_temperature(0)
        liquid_ice.set_heat_latent(335000)
        th = ste2.Stefan(liquid_ice)
        zf = th.monophasic_steady(Tc, z0=2.2e-3)

        plt.plot(
            np.arange(len(zf))*th.dt,  # + popt[1]
            zf,
            '--k',
            label=(
                r'Theoretical front $\delta$ = {:.1e}'.format(
                    np.sqrt(2*-Tc*k_s/(rho_s*Lf))
                )
            )
        )
        """

        plt.grid(True)
        plt.legend(fancybox=True)

        plt.xlabel('Time (s)')
        plt.ylabel('Distance (m)')

        plt.tight_layout()

        plt.figure(figsize=[8, 4.5])
        ax = plt.axes()
        ax.plot(t_space/fps, zf_std/px_mm*1e3, 'sb', mfc='none')
        ax.set_ylabel('t (s)')
        ax.set_ylabel('zf_std (mm)')
        ax.grid(True)

        plt.show()

        # return -Tc-dT_poly

    def define_geometry(self):
        ref_pente = float(self.datas.iloc[0]['alpha'])
        temp = self.current_image_number
        self.read_image(0)
        im = self.current_image.rotate(angle=-ref_pente)
        self.read_image(temp)

        fig, ax = plt.subplots(1, figsize=[10, 7])
        ax.imshow(im, cmap='gray')

        dr = drag.DraggableRectangle(ax)
        plt.show()

        if dr.valid:
            rect = dr.get_rect()
            r0 = int((rect[1] - rect[0])/2)
            rc = int(r0 + rect[0])
            z0 = int(rect[3] - rect[2])

            print('\tz0 : ', z0)
            print('\tr0 : ', r0)
            print('\trc : ', rc)
            print('Please add this value in the table')

    def unfrosting_front(self, zf):
        time_ref = int(self.datas['t_nuc_calc'].values)
        px_mm = 140000  # int(self.datas['px_mm'].values)
        fps = float(self.datas.iloc[0]['fps'])
        try:
            Tc = float(self.datas.iloc[0]['Tc_nuc'])
        except AttributeError:
            Tc = float(self.datas.iloc[0]['Tc_set'])
        ref_pente = float(self.datas.iloc[0]['alpha'])

        z0 = self.datas['z0'].values
        r0 = self.datas['r0'].values
        rc = self.datas['rc'].values

    def get_contour(self, fac=1, plot=True):
        """
        Determine the contour of the current images.

        Please provide x, ymin and ymax from set_contour_ref().
        """
        gradient_im = self.current_image.gradient(
            type='sobel', size=5
        )

        plt.figure()
        plt.imshow(gradient_im, cmap='gray')
        plt.title('gradient image')

        thres = np.median(gradient_im) + 3/4*np.std(gradient_im)
        im_bw = gradient_im
        im_bw[im_bw < thres] = 0
        im_bw[im_bw > 0] = 1

        plt.figure()
        plt.imshow(im_bw, cmap='gray')
        plt.title('binary image')

        im_er = morphology.binary_erosion(
            np.array(im_bw, dtype=bool),
            morphology.square(2)
        )

        plt.figure()
        plt.imshow(im_er, cmap='gray')
        plt.title('first erosion image')

        im_fill = np.array(im_er, dtype=np.int)
        im_fill = morph.binary_fill_holes(im_er)
        # im_fill = np.array(im_fill, dtype=np.int)
        # im_fill = morphology.binary_fill_holes(im_er).as_type(int)

        plt.figure()
        plt.imshow(im_fill, cmap='gray')
        plt.title('fill holes image')

        im_rm = morphology.remove_small_objects(
            im_fill, 6
        )

        plt.figure()
        plt.imshow(im_rm, cmap='gray')
        plt.title('remove small objects < 6')

        im_rm = np.array(im_rm, dtype=np.int16)
        labels = measure.label(im_rm)
        regions = measure.regionprops(labels)

        markers = np.array([r.centroid for r in regions]).astype(np.uint16)
        marker_image = np.zeros_like(im_rm, dtype=np.int64)
        marker_image[markers[:, 0], markers[:, 1]] = np.arange(len(markers)) + 1

        distance_map = ndimage.distance_transform_edt(1 - im_rm)

        filled = morphology.watershed(1 - distance_map, markers=marker_image)

        filled_connected = measure.label(filled != 1, background=0) + 1

        filled_regions = measure.regionprops(filled_connected)
        mean_area = np.mean([r.area for r in filled_regions])
        filled_filtered = filled_connected.copy()
        for r in filled_regions:
            if r.area < 0.25 * mean_area:
                coords = np.array(r.coords).astype(int)
                filled_filtered[coords[:, 0], coords[:, 1]] = 0

        f, (ax0, ax1, ax2) = plt.subplots(1, 3)
        ax0.imshow(im_rm, cmap='gray')
        ax1.imshow(filled_filtered, cmap='Blues')
        ax2.imshow(distance_map, cmap='gray')

        im_er2 = morphology.binary_erosion(
            image=im_fill,
            selem=morphology.square(4)
        )

        plt.figure()
        plt.imshow(im_er2, cmap='gray')
        plt.title('2nd erosion image')

        im_dilate = morphology.binary_dilation(
            im_er2,
            morphology.square(2)
        )

        plt.figure()
        plt.imshow(im_dilate, cmap='gray')
        plt.title('last dilatation image')

        Cx = Cy = np.array([])
        for i in range(0, self.current_image.size[1]-4, 1):
            it = im_dilate[:, i]
            if np.argmax(it) > 0:
                Cx = np.append(Cx, i)
                Cy = np.append(Cy, np.argmax(it))

        if plot:
            fig, ax = plt.subplots(figsize=(11, 7))
            fig.canvas.set_window_title(
                'Contour visualisation - image {:s}'.format(
                    str(self.current_image_number)
                )
            )

            ax.imshow(self.current_image.image, cmap='gray')
            ax.plot(
                Cx,
                Cy,
                '.g', markersize='3',
                label='Contour2 detected'
            )

            plt.legend(fancybox=True, shadow=True)
            plt.grid(True)
            plt.show()

        return Cx, Cy

    def get_contour2(self):
        import cv2
        im = self.current_image.image
        plt.imshow(im, cmap='gray')
        print(np.max(im), np.min(im))

        im = np.array(self.current_image.image*255, dtype=np.uint8)
        cv2.imshow("Image", im)
        cv2.waitKey(0)

        print(np.mean(im))

        ret, thresh = cv2.threshold(im, np.uint8(np.mean(im)), 255, 0)

        cv2.imshow("threshold image", thresh)
        cv2.waitKey(0)

        im2, contours, hierarchy = cv2.findContours(
            thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE
        )

        cnt = contours[0]
        print(cnt)
        cv2.drawContours(im, contours, -1, (0, 255, 0), 3)
        cv2.waitKey(0)

    def get_contours(self, fac=1):
        """Determine contour for all images."""
        temp = self.current_image_number
        widgets = ['Detecting contour',
                   ' ', progressbar.Percentage(),
                   ' ', progressbar.Bar('=', '[', ']'),
                   ' ', progressbar.ETA(),
                   ' ', progressbar.FileTransferSpeed()]
        pbar = progressbar.ProgressBar(
            widgets=widgets, maxval=len(self.image_list)
        )
        pbar.start()

        for i in range(len(self.image_list)):
            self.read_image(i, _hist=False)
            Cx, Cy = self.get_contour(fac=fac, _hist=False)
            if i < 9:
                np.save(
                    self.contour_directory + '/000' + str(i+1) + '_Cx',
                    Cx
                )
                np.save(
                    self.contour_directory + '/000' + str(i+1) + '_Cy',
                    Cy
                )
            elif i < 99:
                np.save(
                    self.contour_directory + '/00' + str(i+1) + '_Cx',
                    Cx
                )
                np.save(
                    self.contour_directory + '/00' + str(i+1) + '_Cy',
                    Cy
                )
            elif i < 999:
                np.save(
                    self.contour_directory + '/0' + str(i+1) + '_Cx',
                    Cx
                )
                np.save(
                    self.contour_directory + '/0' + str(i+1) + '_Cy',
                    Cy
                )
            pbar.update(i)
        pbar.finish()
        self.read_image(temp)

    def manual_contour(self):
        """Detect the contour of the drop over the time."""
        # pente = np.zeros(self.n_image_tot)
        # y0 = np.zeros((self.n_image_tot, self.current_image.size[1]))
        # temp = self.current_image_number

        time_ref = int(self.datas.iloc[0]['t_ref_calc'])
        time_end = int(self.datas.iloc[0]['t_end_calc'])
        if time_end + 4 < self.n_image_tot:
            time_end += 4
        else:
            time_end = self.n_image_tot

        ref_pente = float(self.datas.iloc[0]['alpha'].replace(',', '.'))

        for i in range(time_ref, time_end):
            self.read_image(i)

    def display_contour(self, ni=-1, ref=False):
        """
        Display image with the contour.

        If ni = -1, displays the current_image, else displays image ni
        """
        if ni > -1:
            self.read_image(ni)

        if self.current_image_number < 9:
            Cx = np.load(
                self.data_directory + 'contour/' + '000' +
                str(self.current_image_number+1) + '_Cx' +
                '.npy'
            )
            Cy = np.load(
                self.data_directory + 'contour/' + '000' +
                str(self.current_image_number+1) + '_Cy' +
                '.npy'
            )
        elif self.current_image_number < 99:
            Cx = np.load(
                self.data_directory + 'contour/' + '00' +
                str(self.current_image_number+1) + '_Cx' +
                '.npy'
            )
            Cy = np.load(
                self.data_directory + 'contour/' + '00' +
                str(self.current_image_number+1) + '_Cy' +
                '.npy'
            )
        elif self.current_image_number < 999:
            Cx = np.load(
                self.data_directory + 'contour/' + '0' +
                str(self.current_image_number+1) + '_Cx' +
                '.npy'
            )
            Cy = np.load(
                self.data_directory + 'contour/' + '0' +
                str(self.current_image_number+1) + '_Cy' +
                '.npy'
            )

        fig = plt.figure(figsize=(11, 7))
        fig.canvas.set_window_title(
            'Image with contour n°' + str(self.current_image_number)
        )
        plt.imshow(self.current_image.image, cmap='gray')
        plt.plot(
            Cx, Cy,
            '.', markersize=5, color='#BB0029', alpha=.4,
            label='Contour'
        )
        plt.legend(shadow=True, fancybox=True, loc='best')
        plt.show()
        print(Cy)
        s = sum(Cy)
        print(s)

    def _rot_inv(self, x_rot, y_rot, angle):
        a = angle * np.pi / 180
        x_rot = np.array(x_rot)
        y_rot = np.array(y_rot)
        h_c = int(round(self.current_image.size[0]/2))
        w_c = int(round(self.current_image.size[1]/2))
        x_rot -= w_c
        y_rot -= h_c
        x = x_rot*np.cos(a) - y_rot*np.sin(a)
        y = y_rot*np.cos(a) + x_rot*np.sin(a)
        x += w_c
        y += h_c
        return x, y


class GetFolderName(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        file = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        self.name = file
