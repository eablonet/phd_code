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
import collections

import numpy as np
from scipy import optimize

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

import peakutils as peaks

# local packages
# --------------
from . import Image
from . import read_online_data as rd

from ..geometrybuilder import LineBuilder as lb
from ..geometrybuilder import DraggrableRectangle as drag
from ..geometrybuilder import Geometry as geom

from packages.design import color as ea_color
from packages.library import main as pers_conv

from ..gui import ClaheWindow as cw


# code
# ----
# usefull mehtods
def log(*x):
    """Own print."""
    print(*x, file=sys.stderr)


class Geometry:
    """Class containing the drop geometry."""

    def __init__(self):
        """Init the class."""
        # r -axis data
        # ------
        self.r0 = None
        self.rr = None
        self.rl = None
        self.rc = None

        # z-axis data
        # ------
        self.z0 = None
        self.zf = None
        self.zb = None

    def isEmpty(self):
        """Return True if geometry has been loaded.

        This mean that rr, rl zb or zf are not None.
        """
        if self.zf is None:
            return True
        else:
            return False

    def _get_circle_interpolation(self):
        """Return the circle interpolation."""
        # Get x, yValues
        # --------------
        x = [
            self.rc,  # rc
            self.rl,  # rl
            self.rr,  # rr
        ]
        y = [
            self.z0,  # z0
            self.zb,  # zb
            self.zb  # zb
        ]

        for i in x+y:
            if i is None:
                raise ValueError(
                    'At least one value is missing among :\n' +
                    'rc, rl, rr, z0, zb'
                )

        x = np.array(x)
        y = np.array(y)

        center_estimation = [np.mean(x), np.mean(y)]

        def calc_radius(xc, yc):
            """Calculate the radius of circle."""
            return np.sqrt((x-xc)**2 + (y-yc)**2)

        def minimize(c):
            """Minimize the residual radius."""
            R = calc_radius(*c)
            return R - R.mean()

        center, ier = optimize.leastsq(minimize, center_estimation)
        xc, yc = center
        R = np.mean(calc_radius(xc, yc))

        return int(xc), int(yc), int(R)

    def get_theta(self):
        """Return the initiale contact angle of the drop."""
        xc, yc, R = self._get_circle_interpolation()
        x = np.arange(self.rl, self.rr)
        y = (yc - np.sqrt(R**2 - (x - xc)**2))
        ca_l = np.arctan(y[1] - y[0])*180/np.pi
        ca_r = np.arctan(y[-1] - y[-2])*180/np.pi

        return ca_l, ca_r

    def get_volume_truncated_sphere(self):
        """Return the initiale volume by spherical cap hypothesis."""
        vol = np.pi/6 * self.z0 * (
            3*self.r0**2 + self.z0**2
        )
        return vol

    def get_vol_calculated(self):
        """Return calculated volume by fit circle."""
        xc, yc, R = self._get_circle_interpolation()
        y = np.arange(self.zb, self.z0)
        xl = xc - np.sqrt(R**2 - (y - yc)**2)
        xr = xc + np.sqrt(R**2 - (y - yc)**2)
        vol = np.pi * np.sum(np.power((xr - xl)/2, 2))

        return vol


class Front:
    """Class contains the front information."""

    def __inti__(self):
        """Init the class."""
        self.x = None
        self.y = None

    def isEmpty(self):
        """Return True if front has been loaded."""
        if self.x is None:
            return True
        else:
            return False


class Contour:
    """Class contains the front information."""

    def __init__(self):
        """Init the class."""
        self.xl = None  # xl[t][rl:rc]
        self.yl = None  # "
        self.xr = None  # "
        self.yr = None  # "

        self.time = 0

        self.rc = 0  # max of xl, or min of xr

    def isEmpty(self):
        """Return True if contour has not been loaded."""
        if None in [self.xl, self.yl, self.xr, self.yr]:
            return True
        else:
            return False

    def set_data(self, xl, yl, xr, yr):
        """Set the data."""
        self.xl = xl
        self.yl = yl
        self.xr = xr
        self.yr = yr

        self.time = len(self.xl)
        self.rc = int(np.max(self.xl[0]))

    def get_volume(self):
        """Calculate the volume over time."""
        vol = []  # init vector volume
        plt.figure()
        for t in range(self.time):
            # construct r(y) interpolation
            # -------------
            y = list(
                range(0, int(np.max(
                    np.max(self.yl[0]) - self.yl[t]))+1
                )
            )

            xl = []
            xr = []
            for j in y:
                xl.append(
                    self.xl[t][np.argmin(abs(np.max(self.yl[0]) - self.yl[t] - j))]
                )
                xr.append(
                    self.xr[t][np.argmin(abs(np.max(self.yl[0]) - self.yr[t] - j))]
                )
            xl, xr = np.array(xl), np.array(xr)
            # calcul the volume
            # ---------------
            volume_calc = (
                np.pi/2 * np.sum(np.power(self.rc-xl, 2)) +
                np.pi/2 * np.sum(np.power(xr-self.rc, 2))
            )
            if t % 50 == 0:
                plt.plot(y, (xl-self.rc), '--og')
                plt.plot(y, (xr-self.rc), '--or')

            vol.append(volume_calc)
        plt.show()
        return vol

    def get_height(self, loc=None):
        """Return the height of the drop at a specific radii (loc).

        input
        -----
            loc: int
                Radii of the heigh to obtain. Default is None, and return the
                height on the center.
        """
        z = []
        for i in range(self.time):
            if loc is None:
                z.append(
                    np.max(self.yl[0]) - self.yl[i][-1]
                )

        return z

    def plot(self, dt=10):
        """
        Display the contour over time.

        input
        -----
            dt: int
                Time space between two plots.
        """
        plt.figure(figsize=(16, 9))
        for i in range(0, self.time, dt):
            plt.plot(
                self.xl[i], np.max(self.yl[0]) - self.yl[i],
                ls='--', color='tab:green'
            )
            plt.plot(
                self.xr[i], np.max(self.yl[0]) - self.yr[i],
                ls='--', color='tab:red'
            )
        plt.grid(True)
        plt.xlabel('r (px)')
        plt.ylabel('z_d (t) (px)')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def plot_height(self, loc=None):
        """Display the height's drop at specific location (loc).

        input
        -----
            loc : int
                Radii of the heigh to obtain. Default is None, and return the
                height on the center.
        """
        z = self.get_height(loc)
        plt.plot(range(self.time), z, '--k')
        plt.grid(True)
        plt.xlabel('time (frame)')
        plt.ylabel('z_d (t) (px)')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()


class DataDict(dict):
    """Defin my own dict."""

    def __getitem__(self, key):
        """Overwrite the getitem method."""
        if key in ['t_nuc_mes', 't_end_mes']:
            t = dict.__getitem__(self, key)
            t = t.split(':')
            return int(t[0])*60 + int(t[1]) + int(t[2])/60
        else:
            return dict.__getitem__(self, key)

    def __setitem__(self, key, value):
        """Overwrite the setitem method."""
        """Return t_nuc_mes value."""
        if key in ['t_nuc_mes', 't_end_mes']:
            if type(value) is float and np.isnan(value):
                dict.__setitem__(self, key, None)
            else:
                dict.__setitem__(self, key, value)
        else:
            dict.__setitem__(self, key, value)


class ImageList(list):
    """Define my own list to acces image value."""

    def __getitem__(self, index):
        """Overwrite the getitem, to have acces from 1, not from 0."""
        if index < 0:
            index = range(1, len(self)+1)[index]
        if index not in range(1, len(self)+1):
            raise IndexError('Value out of the range')

        return list.__getitem__(self, index-1)


class Data:
    """Data class. Store all the information about the experiment."""

    def __init__(self):
        """Init the class."""
        self.data = DataDict({
            # identification
            # -------------
            'date': None, 'serie': None,

            # geometry
            # --------
            'r0': None, 'rl': None, 'rc': None, 'rr': None,
            'zb': None, 'z0': None, 'zf': None,
            'ca_l': None, 'ca_r': None,
            'isSymm': None,
            'px_mm': None,

            # time
            # ------
            't_nuc_mes': None, 't_end_mes': None,  # format : hh:mm:ss
            't_nuc_calc': None, 't_ref_calc': None, 't_end_calc': None,
            # format : frame number

            # lieu
            # -----
            'chamber': None,  # bool (1 oui, 0, none)

            # etats
            # -----
            'treated': None,  # bool (1 si traité 0 sinon)
            'infos': None,  # int [0..5]

            # chamber, ambiant and cryostat
            # ----------
            'Ta_cons': None, 'phi_cons': None,
            'Ta_chamber': None, 'phi_chamber': None,
            'Ta_therm': None, 'phi_hydro': None,
            'Tc_sta': None, 'Tc_set': None, 'Tc_nuc': None,
            'ref_therm': None,

            # substrat
            # --------
            'alpha': None, 'ref': None,

            # drop
            # ----
            'n_freeze': None, 'type_freeze': None,
            'insulation': None, 'water_nfo': None,

            # camera
            # ------
            'delay': None, 'exp': None, 'fps': None,
            'delay_exp_real': None, 'fps_real': None,

            # observation
            # --------
            'remarque': None
        })

    def __repr__(self):
        """Print."""
        return 'Experiment id : ' + self.id

    def set_data(self, data):
        """Assign all data values."""
        switch = {
            # id
            # --
            'date': 'date', 'serie': 'serie',

            # geometry
            # --------
            'r0': 'r0', 'rl': 'rl', 'rc': 'rc', 'rr': 'rr',
            'zb': 'zb', 'z0': 'z0', 'zf': 'zf',
            'ca_l': 'ca_l', 'ca_r': 'ca_r',
            'symmetrical': 'isSymm',  # bool : if is symmetrical

            # time
            # ------
            't_nuc_mes': 't_nuc_mes', 't_end_mes': 't_end_mes',
            # format : strint(hh:mm:ss)
            't_nuc_calc': 't_nuc_calc',  # format : frame number
            't_ref_calc': 't_ref_calc',  # format : frame number
            't_end_calc': 't_end_calc',  # format : frame number

            # lieu
            # -----
            'lieu': 'chamber',  # bool (1 oui, 0, none)

            # etats
            # -----
            'completed_2': 'treated',  # bool (1 si traité 0 sinon)
            'etat_2': 'infos',  # int [0..5]

            # chamber, ambiant and cryostat
            # ----------
            'Ta': 'Ta_cons', 'phi': 'phi_cons',
            'Ta_int': 'Ta_chamber', 'phi_int': 'phi_chamber',
            'Ta_ext': 'Ta_therm', 'phi_ext': 'phi_hydro',
            'Tc_sta': 'Tc_sta', 'Tc_set': 'Tc_set',
            'Tc_nuc': 'Tc_nuc',
            'ref_therm': 'ref_therm',

            # substrat
            # --------
            'alpha': 'alpha', 'ref': 'ref',

            # drop
            # ----
            'n_freezing': 'n_freeze', 'type': 'type_freeze',
            'insulation': 'insulation', 'water': 'water_nfo',

            # camera
            # ------
            'delay': 'delay', 'exp': 'exp', 'fps': 'fps',
            'delay_exp_real': 'delay_exp_real', 'real_fps': 'fps_real',

            # observation
            # --------
            'Remarques': 'remarque'
        }
        for s in switch:
            if s in data:
                self.data[switch[s]] = data[s]

    @property
    def id(self):
        """Return unique id of the experiment."""
        if self.data['date'] is None or self.data['serie'] is None:
            return 'None'
        else:
            return self.data['date'] + 'n' + str(self.data['serie'])


class Stack(object):
    """This class allow you to treat image by block.

    Parameters
    ----------
    data_directory : str
        directory containing images

    """

    def __init__(self):
        """Init the Class.

        input
        -----
        date: str
            Choose the date of the manip. Default None don't read any folder.
        serie: int
            Choose the serie of the manip. Default None don't read any folder.

        Notice
        ------
        If date is not None, but serie is None, will read the first serie (1).

        """
        self.path = None
        self.geom = Geometry()

        self.data = Data()

        self.isTreated = False
        self.listResult = []
        self.front = Front()  # x, y data fro front
        self.contour = Contour()  # xl, yl, xr, yr for contour
        self.geometry = Geometry()  # dict with rr, rl, rc, zf, z0 & zb

        self.image = ImageList()  # array containing all Image
        self.range_images = range(0)
        self.n_selected = None

    def __str__(self):
        """Return n# image info."""
        return "Nombre d'image importées : " + str(len(self))

    def __repr__(self):
        """Return n# image info."""
        return self.__str__()

    def __len__(self):
        """Get the number of image."""
        return len(self.image)

    def read_by_date(self, date, serie):
        """Choose the path.

        If the experiement is not in the db, an error occur.

        Parameters
        ----------
        date: str
            select the date of the experiment
        serie : int
            select the serie of the experiment

        """
        df = rd.get_data()  # read data in googlesheet.

        for _, row in df.iterrows():
            cond = (
                row['date'] == date and
                row['serie'] == serie
            )
            if cond:
                self.path = (
                    "/Volumes/EMERYK_HD/0_data/{:s}/n{:#}/".format(
                        date, int(serie)
                    )
                )
                ind = dict(row)
                self.data.set_data(ind)

                # allocate data Values
                # ---------------------
                self.geom.r0 = self.data.data['r0']
                self.geom.rc = self.data.data['rc']
                self.geom.z0 = self.data.data['z0']

        if os.path.exists(self.path + 'result'):
            self.isTreated = True
            self.listResult = [
                l.split('/')[-1].split('.')[0]
                for l in glob.glob(self.path + 'result/*.npy')
            ]
            if 'front' in self.listResult:
                front = np.load(
                    self.path + 'result/front.npy', allow_pickle=True
                )
                self.front.x = front[0]
                self.front.y = front[1]

            if 'contour' in self.listResult:
                log("I pass here")
                contour = np.load(
                    self.path + 'result/contour.npy', allow_pickle=True
                )
                self.contour.set_data(*contour)
                self.contour.plot(20)
                self.contour.plot_height()
                vol = self.contour.get_volume()
                plt.plot(
                    range(len(vol)), [i/vol[0] for i in vol], '--k'
                )

            if 'geom' in self.listResult:
                ge = np.load(
                        self.path + 'result/geom.npy', allow_pickle=True
                    )[()]
                self.geom.rl = ge['rl']
                self.geom.rr = ge['rr']
                self.geom.rc = ge['rc']
                self.geom.zb = ge['zb']
                self.geom.z0 = ge['z0']
                self.geom.zf = ge['zf']
                self.data.data['z0'] = ge['z0']

        if self.path is None:
            raise ValueError(
                "This experiment may not be in the database" +
                ", please add it before"
            )

        self._set_directory()

    def read_select_folder(self):
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

    def _set_directory(self):
        """Define direcoties."""
        self.exporting_directory = self.path + 'export_image/'

    def load_images(self, regex='cam1_*', ext='.tif'):
        """Update the lists of image."""
        # check what exists in data_folder (aka. pretreated images)
        # ----------
        if os.path.isdir(self.path + '_export'):
            listdir = os.listdir(self.path + '_export')

            # ask for what to read
            print('\t' + str(0) + '. ' + 'Original reader image')
            c = 1
            for d in listdir:
                print('\t' + str(c) + '. ' + d)
                c += 1
            num = int(input('Select a folder to read : '))

            # if return is 0, read original image, else read image in the
            # selected folder
            if num > 0:
                image_list = sorted(glob.glob((
                    self.path + '_export/' + listdir[num-1] + '/' +
                    regex + ext
                )))
            else:
                image_list = sorted(glob.glob((
                    self.path + regex + ext
                )))
        # else read directly image
        else:
            # image list
            image_list = sorted(glob.glob((
                self.path + regex + ext
            )))

        self.range_images = range(1, len(image_list)+1)

        widgets = ['Loading images',
                   ' ', progressbar.Percentage(),
                   ' ', progressbar.Bar('=', '[', ']'),
                   ' ', progressbar.ETA(),
                   ' ', progressbar.FileTransferSpeed()]
        pbar = progressbar.ProgressBar(
            widgets=widgets, maxval=len(self.range_images)+1
        )
        pbar.start()
        for n in self.range_images:
            self.image.append(Image.Image())
            self.image[n].read_by_path(image_list[n-1])
            pbar.update(n)
        pbar.finish()
        self.n_selected = 1  # from 1 to range_image

        # check previous treatment
        # ------------
        if os.path.exists(self.path + 'result/crop_data.npy'):
            crop_data = np.load(self.path + 'result/crop_data.npy')
            widgets = ['Is croping all images',
                       ' ', progressbar.Percentage(),
                       ' ', progressbar.Bar('=', '[', ']'),
                       ' ', progressbar.ETA(),
                       ' ', progressbar.FileTransferSpeed()]
            pbar = progressbar.ProgressBar(
                widgets=widgets, maxval=len(self)+1
            )
            pbar.start()
            for i in self.range_images:
                self.image[i].crop(*crop_data)
                pbar.update(i)
            pbar.finish()

    def print_info(self):
        """Show info about the stack."""
        print()
        print('Working directory')
        print('-'*len('Working directory'))
        print(self.path)
        print()
        print('image info')
        print('-'*len('image info'))
        print('Number of image : {0}'.format(len(self)))
        print('Selected selected : {0}'.format(self.n_selected))
        print()
        print('Existing data')
        print('-'*len('Existing data'))
        print('There is results : {0}'.format(self.isTreated))
        if self.isTreated:
            [print('\t{0}'.format(l)) for l in self.listResult]
        print()

    def get_image(self, n=None):
        """Return the expected Image()."""
        if n is None:
            n = self.n_selected

        if n not in self.range_images:
            raise IndexError('{0} is out of the range.'.format(
                n
            ))
        return self.image[n]

    def select_image(self, ni):
        """Change the selected image."""
        if ni in self.range_images:
            self.n_selected = ni
        else:
            raise IndexError('Value out of the range')

    def show_image(self):
        """Show image."""
        self.image[self.n_selected].show_image('Image #{0}/{1}'.format(
            self.n_selected,
            len(self)
        ))

    def crop(self):
            """Crop with DraggableRectangle."""
            # show the image
            # --------------
            plt.figure(figsize=(16, 9))
            plt.gcf().canvas.set_window_title('Crop image')
            ax = plt.subplot(1, 1, 1)
            ax.imshow(self.image[-1].image, cmap='gray')

            # generate rectangle
            # --------------
            dr = drag.DraggableRectangle(ax)
            plt.show()

            # lauch recursive crop
            # --------------
            if dr.valid:
                widgets = ['Is croping all images',
                           ' ', progressbar.Percentage(),
                           ' ', progressbar.Bar('=', '[', ']'),
                           ' ', progressbar.ETA(),
                           ' ', progressbar.FileTransferSpeed()]
                pbar = progressbar.ProgressBar(
                    widgets=widgets, maxval=len(self)+1
                )
                pbar.start()
                for i in self.range_images:
                    self.image[i].crop(*dr.get_rect())
                    pbar.update(i)
                pbar.finish()

                if os.path.exists(self.path + 'result/crop_data.npy'):
                    rect = np.load(self.path + 'result/crop_data.npy')
                    n_rect = dr.get_rect()
                    rect[1] = rect[0] + n_rect[1]
                    rect[0] += n_rect[0]
                    rect[3] = rect[2] + n_rect[3]
                    rect[2] += n_rect[2]
                    np.save(self.path + 'result/crop_data.npy', rect)
                else:
                    np.save(self.path + 'result/crop_data.npy', dr.get_rect())

    def clahe(self):
        """Apply clahe thanks GUI."""
        # get the limit and Size
        # -------
        qapp = QApplication.instance()
        if qapp is None:
            qapp = QApplication(sys.argv)
        app = cw.ClaheWindow(self.image[self.n_selected])
        app.show()
        qapp.exec_()

        # apply to all image (+progressbar)
        # --------
        widgets = ['Applying CLAHE to all images',
                   ' ', progressbar.Percentage(),
                   ' ', progressbar.Bar('=', '[', ']'),
                   ' ', progressbar.ETA(),
                   ' ', progressbar.FileTransferSpeed()]
        pbar = progressbar.ProgressBar(
            widgets=widgets, maxval=len(self)+1
        )
        pbar.start()
        for i in self.range_images:
            self.image[i].equalize_hist_by_clahe(*app.get_data())
            pbar.update(i)
        pbar.finish()

    def rotate(self):
        """Rotate all the images."""
        # apply to all image (+progressbar)
        # --------
        widgets = ['Rotating images',
                   ' ', progressbar.Percentage(),
                   ' ', progressbar.Bar('=', '[', ']'),
                   ' ', progressbar.ETA(),
                   ' ', progressbar.FileTransferSpeed()]
        pbar = progressbar.ProgressBar(
            widgets=widgets, maxval=len(self)+1
        )
        pbar.start()
        for i in self.range_images:
            self.image[i].set_rotate(self.data.data['alpha'])
            pbar.update(i)
        pbar.finish()

    def get_geometry_by_rectangle(self):
        """Get the geometry."""
        alpha = -float(self.data.data['alpha'])
        im = self.image[1].get_rotate(angle=alpha)

        plt.figure(figsize=(16, 9))
        ax = plt.subplot(1, 1, 1)
        ax.imshow(im, cmap='gray')

        dr = drag.DraggableRectangle(ax)
        plt.show()

        if dr.valid:
            rect = dr.get_rect()
            r0 = int((rect[1] - rect[0])/2)
            rc = int(r0 + rect[0])
            z0 = int(rect[3] - rect[2])

            print('Measure')
            print('-'*len('Measure'))
            print('\tz0 : ', z0)
            print('\tr0 : ', r0)
            print('\trc : ', rc)
            print('*Please add this value in the table')

    def get_spatio_col(self, row):
        """Return the diagram spatio temporel at a specific row.

        A spatio is the juxtaposition of column intensity over time.
        We obtain an image of y over t, not y over x.
        """
        sp = np.zeros((self.image[1].get_height(), self.range_images[-1]))
        for n in self.range_images:
            im = self.get_image(n)
            it = im.get_col(row)
            sp[:, n-1] = it
        return sp

    def get_spatio_row(self, col):
        """Return the spatio-temporel diagram at specific column.

        A spatio is the juxtaposition of column intensity over time.
        We obtain an image of y over t, not y over x.

        input
        ----
            col : int
                The column to s
        """
        sp = np.zeros((self.iamge[1].get_height(), self.range_images[-1]))
        for n in self.range_images:
            im = self.get_image(n)
            it, _ = im.col_intensity(col)
            sp[:, n-1] = it
        return sp

    def tracker(self, folder_name='front', pts=[], cmp='pink'):
        """Get the front manually.

        Steps :
        -------
            1. get geometry -- rc, rl, rr, zb, z0  & zf
            2. get volume thanks to double interpolation time + space
        """
        time_ref = int(self.data.data['t_ref_calc'])
        # time to start detection
        time_end = int(self.data.data['t_end_calc'])
        # time to end detection

        # Step 1 - Get geometry
        # ------------
        # generic method to get/save geometry
        def set_geom(geom=None):
            image = [
                self.get_image(time_ref),
                self.get_image(time_end),
            ]
            geom = lb.Geometry(
                image,
                geom=geom,
                px_mm=self.data.data['px_mm'],
            )
            plt.show()
            geom = geom.get_geom()
            np.save(self.path + 'result/geom', geom)
            return geom

        if not self.geom.isEmpty():
            print("There's already geometry data,")
            print("Do you want to overwrite them ? (y/n)")
            a = input()
            if a == 'y':
                geom = set_geom(self.geom)
            elif a == 'n':
                geom = np.load(
                    self.path + 'result/geom.npy', allow_pickle=True
                )[()]
            else:
                raise ValueError("y or n are expected as answer")
        else:
            geom = set_geom()

        self.geom.rr = geom['rr']
        self.geom.rl = geom['rl']
        self.geom.rc = geom['rc']
        self.geom.zf = geom['zf']
        self.geom.z0 = geom['z0']
        self.geom.zb = geom['zb']

        # Step 2 - Double interpolation to detect contour
        # ----------------
        def set_contour(point=None):
            data = {
                't_nuc': time_ref,
                't_end': time_end,
                'geom': geom,
                'range_r': range_r
            }
            image = []
            image.append(self.get_image(time_ref))
            image.append(self.get_image(int((time_ref+time_end)/2)))
            image.append(self.get_image(time_end))
            spatio = lb.SpatioContourTrack(sp, image, data, point=point)
            plt.show()
            contour = spatio.get_contour()  # get the contour over time
            self.contour.set_data(*contour)
            points = spatio.get_tracked_points()

            # save contour
            np.save(self.path + 'result/contour', contour)
            np.save(self.path + 'result/contour_pt', points)

            return contour

        range_r = [
            int(geom['rl'] + 1/3*(geom['rc'] - geom['rl'])),
            int(geom['rl'] + 2/3*(geom['rc'] - geom['rl'])),
            int(geom['rc']),
            int(geom['rc'] + 1/3*(geom['rr'] - geom['rc'])),
            int(geom['rc'] + 2/3*(geom['rr'] - geom['rc'])),
        ]
        sp = {}
        for r in range(len(range_r)):
            sp[r] = self.get_spatio_col(range_r[r])

        if not self.contour.isEmpty():  # il y a un contour
            print("There's already a contour detected,")
            print("Do you want to overwrite them ? (y/n)")
            a = input()
            if a == 'y':  # on réécrit dessus
                contour_pt = np.load(
                    self.path + 'result/contour_pt.npy', allow_pickle=True
                )
                contour = set_contour(
                    point=contour_pt
                )

            elif a == 'n':   # on le garde tel quel
                contour = np.load(
                    self.path + 'result/contour.npy', allow_pickle=True
                )
        else:  # pas encore de contour, on le detect
            contour = set_contour()

        # display the contour
        self.contour.plot()

        # Step 3 - double interpolation to get front
        # ---------------
        frontTracker = lb.SpatioFrontTrack(sp, image, data, contour)
        plt.show()
        front = frontTracker.get_front()

        if not os.path.isfile(self.data_directory + 'result/' + 'front_x.npy'):
            # save front (xfront, yfront)
            np.save(
                self.data_directory + folder_name + '/front',
                front
            )
            # save contour (cx_left, cy_left, cx_right, cy_right)
            np.save(
                self.data_directory + folder_name + '/contour',
                contour
            )
            # save rc, rl, rr, zb, z0, zf
            np.save(
                self.data_directory + folder_name + '/geom',
                geom
            )

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
        x_space = np.arange(len(zf[0, :]))
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

    def get_rfront_loc(self, loc='center'):
        zf = self.read_data()
        t_space, zf, zf_mean, zf_std = self.get_dynamic_front(zf)

        time = np.arange(0, len(zf_mean))
        r = []

        if loc == 'center':
            r0 = self.get_data(int(self.datas['r0'].values))
            for t in time:
                r.append(zf[t, r0])
        elif loc == 'min':
            for t in time:
                r.append(zf[t, 0])
        elif loc == 'max':
            for t in time:
                r.append(zf[t, -1])

        return r


class GetFolderName(QWidget):
    def __init__(self):
        QWidget.__init__(self)
        file = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        self.name = file


if __name__ == '__main__':
    # data example
    # ------------
    data = Data()
    print(data)

    data.data['t_nuc_mes'] = '00:12:30'
    print(data.data['t_nuc_mes'])
    data.data['date'] = '11-07-2019'
    data.data['serie'] = 1
    print(data)
    print(data.data['rr'])

    # stack example
    # -------------
    # instance stack
    stack = Stack()
    print(stack)

    # read an experiment
    stack.read_by_date('17-10-2018', 1)
    stack.print_info()

    # load the images
    stack.load_images()
    stack.print_info()
    stack.select_image(23)
    stack.show_image()

    # make treatement
    stack.rotate()
    stack.crop()
    stack.show_image()
    stack.clahe()
    stack.show_image()

    # pre-info on drop
    stack.get_geometry_by_rectangle()

    # tracking drop geom/contour & front
    stack.tracker()
