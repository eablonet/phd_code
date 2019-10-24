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
IMFT and Toulouse INP

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

from packages.gui import ploter as pl

from packages.design import color as ea_color
from packages.library import main as pers_conv

from ..gui import ClaheWindow as cw

from packages.theory import Stefan as ste


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
        # because zf is only calculated by track_geometry.
        # Is not in the database
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

    def __init__(self):
        """Init the class."""
        self.x = None
        self.y = None

    def isEmpty(self):
        """Return True if front has been loaded."""
        if self.x is None:
            return True
        else:
            return False

    def plot(self, dt=10):
        """Display front plot."""
        plt.figure(figsize=(9, 4.5))
        for i in range(0, len(self.x), dt):
            plt.plot(
                self.x[i], self.y[i],
                ls='--',
            )
        plt.grid(True)
        plt.xlabel('r (px)')
        plt.ylabel('z_d (t) (px)')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def get_x_real(self, rc, pxmm=None):
        """Return true value of x, center on rc.

        input
            rc: int or float
                in px
            pxmm: float or None
                if None return value in px else in mm. Default is None.
        """
        x = []
        for t in range(len(self.x)):
            x.append([])
            for i in range(len(self.x[t])):
                x[t].append(
                    rc - self.x[t][i]
                )
                x[t][i] /= pxmm if pxmm is not None else 1

        return x

    def get_y_real(self, y0=None, pxmm=None):
        """Return true value of y, regard to baseline.

        input
            y0: int or float
                in px. baseline position. If None, the first value of front
                is used.
            pxmm: float or None
                if None return value in px else in mm. Default is None.
        """
        y = []
        if y0 is None:
            c = 0
            while len(self.y[c]) == 0:
                c += 1
            y0 = np.nanmean(self.y[c])

        for t in range(len(self.y)):
            y.append([])
            for i in range(len(self.y[t])):
                y[t].append(y0-self.y[t][i])
                y[t][i] /= pxmm if pxmm is not None else 1

        return y

    def get_front_rloc(self, row):
        """Return the front position at specific location."""
        y = []
        for t in range(len(self.x)):
            for i in range(len(self.x[t])):
                if self.x[t][i] == row:
                    y.append(self.y[t][i])
        return y

    def get_tz0(self):
        None


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
        if any([t is None for t in [self.xl, self.yl, self.xr, self.yr]]):
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
                    self.xl[t][np.argmin(abs(
                        np.max(self.yl[0]) - self.yl[t] - j
                    ))]
                )
                xr.append(
                    self.xr[t][np.argmin(abs(
                        np.max(self.yl[0]) - self.yr[t] - j
                    ))]
                )
            xl, xr = np.array(xl), np.array(xr)
            # calcul the volume
            # ---------------
            volume_calc = (
                np.pi/2 * np.sum(np.power(self.rc-xl, 2)) +
                np.pi/2 * np.sum(np.power(xr-self.rc, 2))
            )

            vol.append(volume_calc)
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

    def plot_volume(self, dt=4, add_to_ax=False, adim=True):
        """Plot the volume over time.

        input
        -----
            add_to_ax : bool
                True to add to an exiting axis, False to plot directly.
                (Display directly the figure)
            adim : bool
                To adimentionalize the volume. Default is true.

        """
        if not add_to_ax:
            fig = pl.Figure()
            fig.set_win_title("Volume over time")
            ax = fig.add_ax(111)
            ax.ax.set_title("Evolution du volume au cours du temps")
            ax.xlabel('Time (frame)')
            ax.ylabel(r'$V/V_0$')

        else:
            ax = plt.gca()

        vol = self.get_volume()

        if adim:
            vol = [v / vol[0] for v in vol]

        ax.plot(
            range(len(vol)), vol,
            ls='-.', marker='o',
            color='tab:red',
            markevery=dt,  mfc='none', mec='tab:red',
        )

        if not add_to_ax:
            fig.show()


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
            'px_mm': 'px_mm',

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

    @property
    def fps(self):
        """Return the fps of the experiment."""
        if self.data['fps_real'] is None:
            return self.data['fps']
        else:
            return self.data['fps_real']

    @property
    def Tnuc(self):
        """Return Nucleation temperature."""
        if self.data['Tc_nuc'] is None:
            return self.data['Tc_set']
        else:
            return self.data['Tc_nuc']


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

        self.stefan = ste.Stefan()
        self.__initStefan__()

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

        if os.path.exists(self.path + 'front'):
            # reader for old data
            # =================
            zf = np.load(self.path + 'front' + '/interpolation.npy')

            x, y = [], []
            for t in range(zf.shape[0]):
                x.append([])
                y.append([])
                for r in range(zf.shape[1]):
                    if zf[t, r] != 0:
                        x[t].append(r)
                        y[t].append(zf[t, r])
            self.front.x = x
            self.front.y = y

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
                contour = np.load(
                    self.path + 'result/contour.npy', allow_pickle=True
                )
                self.contour.set_data(*contour)

            if 'geom' in self.listResult:
                ge = np.load(
                        self.path + 'result/geom.npy', allow_pickle=True
                    )[()]

                # set geom data
                # -------------
                self.geom.rl = ge['rl']
                self.geom.rr = ge['rr']
                self.geom.rc = ge['rc']
                self.geom.zb = ge['zb']
                self.geom.z0 = ge['z0']
                self.geom.zf = ge['zf']
                self.data.data['z0'] = ge['z0']

                # update Stefan values
                # -------------
                self.stefan.set_geometry(
                    H=[
                        0,
                        (self.geom.zb - self.geom.z0)/self.data.data['px_mm']*1e3,
                        0
                    ],
                    dz=[0, .005, 0],
                    unit='mm'
                )
        else:
            os.mkdir(self.path + 'result')

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

    def __initStefan__(self):
        """Init Stefan solver."""
        self.stefan.solver = [0, 0, 2, 0]
        self.stefan.dilatation = True
        self.stefan.boundaries = [0, 0]  # bottom, top
        self.stefan.set_geometry(dz=[0, .05, 0], unit='mm')
        self.stefan.set_time_condition(ts=0, tend=100, dt=.01, auto=True)
        # times

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
            self.image[i].set_rotate(-self.data.data['alpha'])
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

    def track_geometry(self):
        """Detect the geomtry."""
        t_ref = int(self.data.data['t_ref_calc'])
        t_end = int(self.data.data['t_end_calc'])

        # Step 1 - Get geometry
        # ------------
        # generic method to get/save geometry
        def set_geom(geom=None):
            image = [
                self.get_image(t_ref),
                self.get_image(t_end),
            ]
            geom = lb.Geometry(
                image,
                geom=geom,
                px_mm=self.data.data['px_mm'],
            )
            plt.show()
            geom = geom.get_geom()
            np.save(self.path + 'result/geom', geom)

            # update geometry data
            # ====================
            self.geom.rr = geom['rr']
            self.geom.rl = geom['rl']
            self.geom.rc = geom['rc']
            self.geom.zf = geom['zf']
            self.geom.z0 = geom['z0']
            self.geom.zb = geom['zb']

            # update Stefan values
            # -------------
            self.stefan.set_geometry(
                H=[0, (self.zb - self.z0)/self.data.data['px_mm'], 0],
                dz=[0, .005, 0],
                unit='mm'
            )

        if not self.geom.isEmpty():
            print("There's already geometry data,")
            print("Do you want to overwrite them ? (y/n)")
            a = input()

            while True:
                if a == 'y':
                    set_geom(self.geom)
                    break
                elif a == 'n':
                    break
                else:
                    print("y or n are expected as answer")
                    a = input()
        else:
            set_geom()

    def track_contour(self):
        """Detect the contour from spatio."""
        if self.geom.isEmpty():
            raise ImportError(
                "Geometry is not available." +
                "Please to track the geomtry first."
            )

        def set_contour(point=None):
            # spatio location
            # ----------
            range_r = [
                int(self.geom.rl + 1/3*(self.geom.rc - self.geom.rl)),
                int(self.geom.rl + 2/3*(self.geom.rc - self.geom.rl)),
                int(self.geom.rc),
                int(self.geom.rc + 1/3*(self.geom.rr - self.geom.rc)),
                int(self.geom.rc + 2/3*(self.geom.rr - self.geom.rc)),
            ]

            # data dictionary
            # -------
            t_ref = int(self.data.data['t_ref_calc'])
            t_end = int(self.data.data['t_end_calc'])
            data = {
                't_nuc': t_ref,
                't_end': t_end,
                'geom': self.geom,
                'range_r': range_r
            }

            # select images
            # ----------
            image = []
            image.append(self.get_image(t_ref))
            image.append(self.get_image(int((t_ref+t_end)/2)))
            image.append(self.get_image(t_end))

            # select spatio
            # ----------
            sp = {}
            for r in range(len(range_r)):
                sp[r] = self.get_spatio_col(range_r[r])
            spatio = lb.SpatioContourTrack(sp, image, data, point=point)
            plt.show()
            contour = spatio.get_contour()  # get the contour over time
            self.contour.set_data(*contour)
            points = spatio.get_tracked_points()

            # save contour
            np.save(self.path + 'result/contour', contour)
            np.save(self.path + 'result/contour_pt', points)

        if not self.contour.isEmpty():  # il y a un contour
            print("There's already a contour detected,")
            print("Do you want to overwrite them ? (y/n)")
            a = input()
            while True:
                if a == 'y':  # on réécrit dessus
                    contour_pt = np.load(
                        self.path + 'result/contour_pt.npy', allow_pickle=True
                    )
                    set_contour(
                        point=contour_pt
                    )
                    break
                elif a == 'n':   # on le garde tel quel
                    break
                else:
                    print("y or n are expected as answer")
                    a = input()

        else:  # pas encore de contour, on le detect
            set_contour()

    def track_front(self):
        """Get the front manually."""
        if self.contour.isEmpty():
            raise ImportError(
                "Geometry is not available." +
                "Please to track the geomtry first."
            )

        def set_front(point=None):
            # create data field to pass to UI
            # ----------
            # spatio location
            # ----------
            range_r = [
                int(self.geom.rl + 1/3*(self.geom.rc - self.geom.rl)),
                int(self.geom.rl + 2/3*(self.geom.rc - self.geom.rl)),
                int(self.geom.rc),
                int(self.geom.rc + 1/3*(self.geom.rr - self.geom.rc)),
                int(self.geom.rc + 2/3*(self.geom.rr - self.geom.rc)),
            ]

            t_ref = int(self.data.data['t_ref_calc'])
            t_end = int(self.data.data['t_end_calc'])

            data = {
                't_nuc': t_ref,
                't_end': t_end,
                'geom': self.geom,
                'range_r': range_r
            }

            # create image and spatio
            # ---------
            image = []
            image.append(self.get_image(t_ref+10))
            image.append(self.get_image(int((t_ref+t_end)/2)))
            image.append(self.get_image(t_end-10))

            sp = {}
            for r in range(len(range_r)):
                sp[r] = self.get_spatio_col(range_r[r])

            # track front
            # -----------
            frontTracker = lb.SpatioFrontTrack(
                sp, image, data, self.contour, point
            )
            plt.show()

            # load data and save them
            # -------------
            front = frontTracker.get_front()
            pt = frontTracker.get_tracked_points()
            np.save(self.path + 'result/front', front)
            np.save(self.path + 'result/front_pt', pt)

            # set self.front points
            # ---------
            self.front.x = front[0]
            self.front.y = front[1]


        if not self.front.isEmpty():
            print("There's already a front detected,")
            print("Do you want to overwrite it ? (y/n)")
            a = input()

            while True:
                if a == 'y':  # on réécrit dessus
                    front_pt = np.load(
                        self.path + 'result/front_pt.npy', allow_pickle=True
                    )
                    set_front(
                        point=front_pt
                    )
                    break
                elif a == 'n':   # on le garde tel quel
                    break
                else:
                    print("y or n are expected as answer")
                    a = input()
        else:
            set_front()

    def front_profile(self, dt=1):
        """View front profile in usi."""
        pxmm = self.data.data['px_mm']
        x, y = np.zeros(len(self.front.x)), np.zeros(len(self.front.x))
        x = self.front.x
        y = self.front.y

        # centering data
        xf = []
        if self.geom.rc is not None:
            for t in range(len(x)):
                xf.append([])
                for i in range(len(x[t])):
                    xf[t].append((self.geom.rc - x[t][i]) / pxmm)
        else:
            for t in range(len(x)):
                xf.append([])
                for i in range(len(x[t])):
                    xf[t].append(x[t][i] / pxmm)

        # remove baseline height
        yf = []
        if self.geom.zb is None:
            c = 0
            while len(y[c]) == 0:
                c += 1
            for t in range(len(y)):
                yf.append([])
                for i in range(len(y[t])):
                    yf[t].append((np.nanmean(y[c])-y[t][i]) / pxmm)
        else:
            for t in range(len(y)):
                yf.append([])
                for i in range(len(y[t])):
                    yf[t].append((self.geom.zb - y[t][i]) / pxmm)
        cmap = plt.get_cmap('coolwarm', len(x))
        fig = pl.Figure()
        ax = fig.add_ax(111)
        for t in range(0, len(x), dt):
            ax.plot(
                xf[t], yf[t],
                ls='--', color=cmap(t)
            )
        time = np.arange(len(x)) / self.data.fps
        fig.colorbar(mapname='coolwarm', values=time, clabel='Time')
        ax.xlabel(r'$r$ (mm)')
        ax.ylabel(r'$z_d$ (t) (mm)')
        # plt.gca().set_aspect('equal', adjustable='box')
        fig.show()

    def get_dynamic_mean_front(self):
        """Return mean dynamic front."""
        pxmm = self.data.data['px_mm']
        if self.geom.zb is None:
            yf = self.front.get_y_real(y0=None, pxmm=pxmm)
        else:
            yf = self.front.get_y_real(y0=self.geom.zb, pxmm=pxmm)

        time = np.arange(0, len(yf))

        y_mean = []
        y_std = []

        for t in time:
            y_mean.append(np.mean(yf[t]))
            y_std.append(np.std(yf[t]))

        return y_mean, y_std

    def get_dynamic_loc_front(self, col):
        """Return dynamic front at a specific radius(row in px)."""
        pxmm = self.data.data['px_mm']
        xf = self.front.x
        if self.geom.zb is None:
            yf = self.front.get_y_real(y0=None, pxmm=pxmm)
        else:
            yf = self.front.get_y_real(y0=self.geom.zb, pxmm=pxmm)

        yi = []
        time = []
        for i, x in enumerate(xf):
            id = np.where(x == col)
            if len(id[0]) == 1:
                time.append(i)
                yi.append(yf[i][id[0][0]])
        return yi, time

    def plot_dynamic_mean_front(self):
        """Display the mean front of solidification."""
        y_m, _ = self.get_dynamic_mean_front()

        dt0 = int(
                self.data.data['t_nuc_calc'] - self.data.data['t_ref_calc'] - 1
            ) / self.data.fps
        time = np.arange(len(y_m)) / self.data.fps

        fig = pl.Figure()
        fig.set_win_title("Mean front dynamic")
        ax = fig.add_ax(111)
        ax.plot(
            time-dt0, y_m,
            marker='.', ls='-.', color='tab:red', markevery=2
        )
        ax.ax.set_title("Evolution du front moyen de solidification")
        ax.xlabel('time (s)')
        ax.ylabel(r'$z_i(m)$')
        fig.show()

    def get_qa(self, init_val=-7000, eps=.1, N=100):
        """Return the closer qa to fit tz0.

        input
        -----
            init_val : float
                Init air-liquid flux value.
            eps : float
                Precision of solving Dtz0 < eps
            N : int
                Maximum step calculation
        """
        tz0 = self.get_tz0(usi=True)

        s = ste.Stefan()
        qa = init_val
        c = 0  # compteur pour itération maximale < N

        # init solver
        # -----------
        s.solver = [0, 1, 1, 0]  # air, liq, ice, subs
        s.dilatation = False
        s.boundaries = [0, 2]  # bottom, top
        s.boundValues = [self.data.Tnuc, qa]  # bottom, top
        z0 = (self.geom.zb - self.geom.z0)/self.data.data['px_mm']
        s.set_geometry(
            H=[0, z0*1e3, 0], dz=[0, .005, 0], unit='mm'
        )  # air, liq, sub
        s.set_time_condition(
            ts=0, tend=self.get_tf(usi=True),
            dt=.1,
            auto=False
        )  # times
        s.solve()

        b0 = 2e12  # ecart of tz0/tz0_ste from 1, init at infinity
        print('Solving qa, please wait...')

        # loop optimize qa
        # ----------------
        while abs(s.tz0 - tz0) > eps and c < N:
            # if stefan.tz0 > tz0 => qa might decrease
            # else stefan.tz0 > tz0 => might increase
            if abs(tz0/s.tz0-1) < abs(b0):
                b1 = tz0/s.tz0-1
            else:
                b1 = b0

            if b1*b0 < 0:
                # changement de signe (on passe d'une surestimation
                # à une sous estimation ou le contraire)
                b0 = b1
                qa1 = qa*(b0+1)
                qa = (qa + qa1)/2
            else:
                # on continue de sous estimer ou du surestimer
                b0 = b1
                qa *= b0+1

            s.boundValues[1] = qa
            s.solve()

            log('\tqa : ', qa)
            log('\ttz0 / tz0.ste : ', tz0, s.tz0)
            c += 1

        return qa

    def get_T_for_airflux(self):
        """Return heat field for flux condition at air-liq interface."""
        # load values
        # -----------
        qa = self.get_qa()
        z0 = (self.geom.zb - self.geom.z0)/self.data.data['px_mm']

        # init solver
        # ----------
        s = ste.Stefan()
        s.solver = [0, 1, 1, 0]  # air, liq, ice, subs
        s.dilatation = False
        s.boundaries = [0, 2]  # bottom, top
        s.boundValues = [self.data.Tnuc, qa]  # bottom, top
        s.set_geometry(
            H=[0, z0*1e3, 0], dz=[0, .005, 0], unit='mm'
        )  # air, liq, sub
        s.set_time_condition(
            ts=0, tend=self.get_tf(usi=True),
            dt=.1,
            auto=False
        )  # times

        # solve
        # -----
        s.solve()

        return s.fields.T

    def get_zi_for_airflux(self):
        """Return ice front position for flux condition at air-liq interface."""
        # load values
        # -----------
        qa = self.get_qa()
        z0 = (self.geom.zb - self.geom.z0)/self.data.data['px_mm']

        # init solver
        # ----------
        s = ste.Stefan()
        s.solver = [0, 1, 1, 0]  # air, liq, ice, subs
        s.dilatation = False
        s.boundaries = [0, 2]  # bottom, top
        s.boundValues = [self.data.Tnuc, qa]  # bottom, top
        s.set_geometry(
            H=[0, z0*1e3, 0], dz=[0, .005, 0], unit='mm'
        )  # air, liq, sub
        s.set_time_condition(
            ts=0, tend=self.get_tf(usi=True),
            dt=.1,
            auto=False
        )  # times

        s.solve()

        return s.fields.T

    def get_tz0(self, usi=False):
        """Return tz0 time."""
        z0 = (self.geom.zb - self.geom.z0) / self.data.data['px_mm']
        yi, _ = self.get_dynamic_mean_front()
        t = range(len(yi))
        tz0 = t[np.argmin(np.abs(z0 - np.array(yi)))]
        if usi:
            tz0 /= self.data.fps
        return tz0

    def get_tz0ste(self):
        """Return theoretical Stefan time."""
        return self.stefan.tz0

    def get_t12(self, usi=False):
        """Return tz0 time."""
        z0 = (self.geom.zb - self.geom.z0) / self.data.data['px_mm']
        yi, _ = self.get_dynamic_mean_front()
        t = list(range(len(yi)))
        t12 = t[np.argmin(np.abs(z0/2 - np.array(yi)))]
        if usi:
            if self.data.data['fps_real'] is None:
                fps = self.data.data['fps']
            else:
                fps = self.data.data['fps_real']
            t12 /= fps
        return t12

    def get_tf(self, usi=False):
        """Return tz0 time."""
        yi = self.front.y
        tf = len(yi)-1
        if usi:
            if self.data.data['fps_real'] is None:
                fps = self.data.data['fps']
            else:
                fps = self.data.data['fps_real']
            tf /= fps
        return tf

    def get_ice_volume(self):
        """Calculate the volume of ice over time."""
        vol = []  # init vector volume

        # get mean front position in px
        yi = self.front.y
        zi = []
        for y in yi:
            zi.append(int(np.mean(y)))

        for t in range(len(zi)):
            # construct r(y) interpolation
            # -------------
            y = list(
                range(0, self.geom.zb - zi[t])
            )

            xl = []
            xr = []
            for j in y:
                xl.append(
                    self.contour.xl[t][np.argmin(abs(
                        self.geom.zb - self.contour.yl[t] - j
                    ))]
                )
                xr.append(
                    self.contour.xr[t][np.argmin(abs(
                        self.geom.zb - self.contour.yr[t] - j
                    ))]
                )
            xl, xr = np.array(xl), np.array(xr)

            # calcul the volume
            # ---------------
            volume_calc = (
                np.pi/2 * np.sum(np.power(self.geom.rc-xl, 2)) +
                np.pi/2 * np.sum(np.power(xr-self.geom.rc, 2))
            )

            vol.append(volume_calc)

        return vol

    def plot_vol_mean_front_comparison(self):
        """Plot the volume from contour and volume from front."""
        vol_d = self.contour.get_volume()
        dt0 = int(
                self.data.data['t_nuc_calc'] - self.data.data['t_ref_calc'] - 1
            )
        v0 = vol_d[dt0]
        dt0 /= self.data.fps
        vol_d = [i / self.data.data['px_mm']**3 for i in vol_d]
        vol_i = self.get_ice_volume()
        vol_tot = [
            i*(1 - 916/999)/self.data.data['px_mm']**3 + v0/self.data.data['px_mm']**3
            for i in vol_i
        ]
        time = np.arange(len(vol_d)) / self.data.fps
        fig = pl.Figure()
        ax = fig.add_ax(111)
        ax.plot(
            time-dt0, vol_d, color='tab:blue',
            markevery=2, marker='s',
            ls='--',
            label='Volume de la goutte (mm^3)'
        )
        ax.plot(
            time-dt0, vol_tot, color='tab:red',
            markevery=2, marker='s',
            ls='--',
            label='Construction du volume total (mm^3)'
        )
        ax.legend()

        fig.show()

    def plot_mean_vs_loc_kinetics(self):
        """Display the mean front and 5 locals front of solidification."""
        fig = pl.Figure()
        fig.set_win_title("Mean front dynamic")
        ax = fig.add_ax(111)
        y_m, _ = self.get_dynamic_mean_front()
        dt0 = int(
                self.data.data['t_nuc_calc'] - self.data.data['t_ref_calc'] - 1
            ) / self.data.fps

        time = np.arange(len(y_m)) / self.data.fps
        ax.plot(
            time-dt0, y_m,
            marker='.', ls='-', color='tab:red', markevery=2,
            label='Mean front'
        )
        range_r = [
            int(self.geom.rl + 1/3*(self.geom.rc - self.geom.rl)),
            int(self.geom.rl + 2/3*(self.geom.rc - self.geom.rl)),
            int(self.geom.rc),
            int(self.geom.rc + 1/3*(self.geom.rr - self.geom.rc)),
            int(self.geom.rc + 2/3*(self.geom.rr - self.geom.rc)),
        ]
        cmap = plt.cm.get_cmap('tab20', 5)
        for i, r in enumerate(range_r):
            yi, time = self.get_dynamic_loc_front(r)
            time = np.array(time)/self.data.fps

            ax.plot(
                time-dt0, yi,
                marker='.', ls='-.', color=cmap(i), markevery=2,
                label='Front location : {:.1f}mm'.format(
                    (r - self.geom.rc)/self.data.data['px_mm']*1e3
                )
            )

        ax.ax.set_title("Mean kinetic, and locals kinetics")
        ax.xlabel('time (s)')
        ax.ylabel(r'$z_i(m)$')
        ax.legend()
        fig.show()

    def fit_mean_kinetic(self, win=.8):
        """Return the best fit.

        There is an automatic windowing to conserve only point > 0 and 80% of
        first remaining points (to elimate last acceleratio due to dilatation)

        input:
        ------
            win: float [0..1]
                Percent of firsts points > 0 to keep in fitting data.
                Default .8
        """
        def fit_t12(x, a):
            """Fit with impose 1/2 slope."""
            return a + x/2

        def fit_best(x, a, b):
            """Fit with impose best exposant."""
            return a*x + b

        # data
        # ----
        y_mean, _ = self.get_dynamic_mean_front()
        time = np.arange(len(y_mean)) / self.data.fps
        dt0 = int(
                1 + self.data.data['t_nuc_calc'] - self.data.data['t_ref_calc']
            ) / self.data.fps

        # log10 conversion
        # ----------------
        idx = np.where(time - dt0 > 0)
        ylog = np.log10(np.array(y_mean)[idx[0]])

        # windowing 0 to 80%
        # -----------------
        tlog = time[idx[0]]-dt0
        idx = np.where(tlog < win*np.max(tlog))
        tlog = np.log10(tlog[idx[0]])
        ylog = ylog[idx[0]]

        # fits
        p12, _ = optimize.curve_fit(
            fit_t12, tlog, ylog
        )
        pbest, _ = optimize.curve_fit(
            fit_best, tlog, ylog
        )
        time_fit = np.linspace(0, np.max(time-dt0), 100)
        y_t12 = 10**(p12[0])*time_fit**(1/2)
        y_best = 10**(pbest[1])*time_fit**(pbest[0])

        # log('fit : ', popt[0])
        # log('10**fit : ', 10**popt[0])
        # log('time_fit : ', time_fit)
        # log('time_fit**(1/2): ', time_fit**(1/2))
        # log('10**(fit)*time_fit**(1/2): ', 10**popt[0]*time_fit**(1/2))

        # display figure
        # --------------
        fig = pl.Figure()
        ax = fig.add_ax(111)
        ax.loglog(
            time-dt0, y_mean,
            '-.o',
            color='tab:red', mfc='none', label='data'
        )
        ax.loglog(
            time_fit, y_t12,
            '-b', color='tab:blue',
            label=r'Fit $\sim t^{1/2}$'
        )
        ax.loglog(
            time_fit, y_best,
            '-', color='tab:green',
            label=r'Best Fit $\sim t^{{{:.2f}}}$'.format(pbest[0])
        )
        ax.legend()
        ax.xlabel('Time (s)')
        ax.ylabel('Height (m)')
        fig.show()

        return p12, pbest
    # to check if to keep or to erase

    def make_montage(self, front=True, contour=True):
        """Create a montage.

        Parameters
        ----------
        front : bool
            Enable front diplaying
        contour : bool
            Enable contour diplaying

        """
        tnuc = self.data.data['t_ref_calc']
        tend = self.data.data['t_end_calc']
        dt = int((tend-tnuc)/5)

        fig = pl.Figure()
        ax = []
        time = np.int64(np.append(tnuc + np.arange(5)*dt, tend))
        cmap = plt.get_cmap('coolwarm', 6)
        for i, t in enumerate(time):
            ax.append(fig.add_ax(231+i))
            ax[-1].imshow(self.get_image(t).image)

            n = int(t - tnuc)
            if front:
                x, y = self.front.x[n], self.front.y[n]
                ax[-1].plot(
                    x, y,
                    ls='--', color=cmap(i), alpha=.5
                )
            if contour:
                x1, y1 = self.contour.xl[n], self.contour.yl[n]
                x2, y2 = self.contour.xr[n], self.contour.yr[n]
                ax[-1].plot(
                    x1, y1, x2, y2,
                    ls='--', color=cmap(i), alpha=.5
                )
            ax[-1].ax.axis('off')
            ax[-1].title('t : {:.1f}s'.format((t - tnuc)/ self.data.fps))
            ax[-1].ax.grid(False)
        fig.show()

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
    # data = Data()
    # print(data)
    #
    # data.data['t_nuc_mes'] = '00:12:30'
    # print(data.data['t_nuc_mes'])
    # data.data['date'] = '11-07-2019'
    # data.data['serie'] = 1
    # print(data)
    # print(data.data['rr'])

    # stack example
    # -------------
    stack = Stack()

    # read an experiment
    # ==========
    stack.read_by_date('22-10-2018', 7)
    # stack.print_info()

    # load the images
    # ==========
    stack.load_images()
    # stack.print_info()
    # stack.select_image(23)
    # stack.show_image()

    # make treatement
    # =========
    # stack.crop()
    # stack.rotate()
    # stack.clahe()
    # stack.show_image()

    # pre-info on drop
    # =========
    # stack.get_geometry_by_rectangle()

    # tracking drop geom/contour & front
    # ========
    # stack.track_geometry()
    # stack.track_contour()
    # stack.track_front()

    stack.contour.plot()
    stack.contour.plot_volume()
    stack.front.plot(dt=1)
    stack.plot_vol_mean_front_comparison()
    stack.front_profile(dt=2)
    stack.plot_dynamic_mean_front()
    stack.plot_mean_vs_loc_kinetics()
    stack.fit_mean_kinetic(.7)
    stack.make_montage()
    log(stack.get_t12(True))
    log(stack.get_tf(True))
    log(stack.get_tz0(True))
    log(stack.get_qa())
