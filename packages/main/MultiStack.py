# multiple stack plotter

# global packages
import numpy as np
import matplotlib.pyplot as plt

# local packages
# -------------
import packages.main.Stack2 as st
import packages.theory.Stefan as ste
import packages.gui.ploter as pl


class Load:
    def __init__(self):
        """Init."""
        self.n =  0
        self.paths = []
        self.stack = []
        self.stefan = []

    def add_path(self, date, serie):
        """Add stacks to db.

        input
        -----
            date: str
                Date de la manip
            serie: int
                numero de la serie
        """
        self.stack.append(st.Stack())
        self.stack[-1].read_by_date(date, serie)

    def kinetic(self, colorby='None'):
        """Plot kinetics."""
        fig = pl.Figure()
        ax = fig.add_ax(111)
        cmap = plt.cm.get_cmap('tab20', len(self.stack))

        for i, s in enumerate(self.stack):
            ym, ystd = s.get_dynamic_mean_front()
            t = np.arange(len(ym)) / s.data.fps
            dt0 = int(
                    s.data.data['t_nuc_calc'] - s.data.data['t_ref_calc']
                ) / s.data.fps

            ax.plot(
                (t-dt0)/s.get_tz0ste(), np.array(ym)/s.get_z0(usi=True),
                c=cmap(i),
                # c=cmap(phi[s.data.data['phi_cons']]),
                ls='none', marker='o', mfc='none',
                # label='__legend__' if leg[s.data.data['phi_cons']] else
                # r'phi = {:.1f}'.format(s.data.data['phi_cons'])
                label='T_w ; {:.0f}'.format(s.data.Tnuc)
            )
            # leg[s.data.data['phi_cons']] = True

        ax.xlabel(r'$t/t_{ste}$')
        ax.ylabel(r'$z/z_0$')
        ax.legend()
        fig.show()

    def kinetic_adim_tz0(self):
        """Plot kinetics."""
        fig = pl.Figure()
        ax = fig.add_ax(111)
        cmap = plt.cm.get_cmap('tab20c', 4)
        phi = {40: 0, 50: 1, 60: 2, 70: 3}
        leg = {40: False, 50: False, 60: False, 70: False}
        for i, s in enumerate(self.stack):
            ym, ystd = s.get_dynamic_mean_front()
            t = np.arange(len(ym)) / s.data.fps
            dt0 = int(
                    s.data.data['t_nuc_calc'] - s.data.data['t_ref_calc']-1
                ) / s.data.fps
            ax.plot(
                (t-dt0)/s.get_tz0(usi=True), np.array(ym)/s.get_z0(),
                c=cmap(phi[s.data.data['phi_cons']]),
                ls='-.', marker='v', mfc='none',
                label='__legend__' if leg[s.data.data['phi_cons']] else
                r'phi = {:.1f}'.format(s.data.data['phi_cons'])
            )
            leg[s.data.data['phi_cons']] = True
        ax.xlabel(r'$t/t_{z_0}$')
        ax.ylabel(r'$z/z_0$')
        ax.legend()
        fig.show()

    def kinetic_adim_tste(self):
        """Plot kinetics."""
        fig = pl.Figure()
        ax = fig.add_ax(111)
        cmap = plt.cm.get_cmap('tab20c', 4)
        phi = {40: 0, 50: 1, 60: 2, 70: 3}
        leg = {40: False, 50: False, 60: False, 70: False}
        for i, s in enumerate(self.stack):
            ym, ystd = s.get_dynamic_mean_front()
            t = np.arange(len(ym)) / s.data.fps
            dt0 = int(
                    s.data.data['t_nuc_calc'] - s.data.data['t_ref_calc']-1
                ) / s.data.fps
            # adim
            # ax.plot(
            #     (t-dt0)/s.get_tz0(), (np.array(ym)/(s.geom.zb - s.geom.z0)*s.data.data['px_mm']),
            #     c=cmap(phi[s.data.data['phi_cons']]),
            #     ls='-.', marker='v', mfc='none',
            #     label='__legend__' if leg[s.data.data['phi_cons']] else
            #     r'phi = {:.1f}'.format(s.data.data['phi_cons'])
            # )
            # usi
            z0 = (s.geom.zb - s.geom.z0)/s.data.data['px_mm']
            ax.plot(
                (t-dt0)/s.get_tz0(usi=True), np.array(ym)/z0,
                c=cmap(phi[s.data.data['phi_cons']]),
                ls='-.', marker='v', mfc='none',
                label='__legend__' if leg[s.data.data['phi_cons']] else
                r'phi = {:.1f}'.format(s.data.data['phi_cons'])
            )
            leg[s.data.data['phi_cons']] = True
        ax.xlabel(r'$t/t_{z_0}$')
        ax.ylabel(r'$z/z_0$')
        ax.legend()
        fig.show()

    def kinetic_adim_tf(self):
        """Plot kinetics."""
        fig = pl.Figure()
        ax = fig.add_ax(111)
        cmap = plt.cm.get_cmap('tab20c', 4)
        phi = {40: 0, 50: 1, 60: 2, 70: 3}
        leg = {40: False, 50: False, 60: False, 70: False}
        for i, s in enumerate(self.stack):
            ym, ystd = s.get_dynamic_mean_front()
            t = np.arange(len(ym)) / s.data.fps
            dt0 = int(
                    s.data.data['t_nuc_calc'] - s.data.data['t_ref_calc']-1
                ) / s.data.fps
            # adim
            # ax.plot(
            #     (t-dt0)/s.get_tz0(), (np.array(ym)/(s.geom.zb - s.geom.z0)*s.data.data['px_mm']),
            #     c=cmap(phi[s.data.data['phi_cons']]),
            #     ls='-.', marker='v', mfc='none',
            #     label='__legend__' if leg[s.data.data['phi_cons']] else
            #     r'phi = {:.1f}'.format(s.data.data['phi_cons'])
            # )
            # usi
            z0 = (s.geom.zb - s.geom.z0)/s.data.data['px_mm']
            ax.plot(
                (t-dt0)/s.get_tz0(usi=True), np.array(ym)/z0,
                c=cmap(phi[s.data.data['phi_cons']]),
                ls='-.', marker='v', mfc='none',
                label='__legend__' if leg[s.data.data['phi_cons']] else
                r'phi = {:.1f}'.format(s.data.data['phi_cons'])
            )
            leg[s.data.data['phi_cons']] = True
        ax.xlabel(r'$t/t_{z_0}$')
        ax.ylabel(r'$z/z_0$')
        ax.legend()
        fig.show()

    def tip_angle(self):
        """Plot tip angle."""
        fig = pl.Figure()
        ax = fig.add_ax(111)
        for s in self.stack:
            alpha = s.contour.tip_angle()
            amin = np.min(alpha)
            ax.plot(s.data.data['phi_cons'], amin, 's')
        fig.show()

    def tz0_repr(self, param='phi_int'):
        fig = pl.Figure()
        ax = fig.add_ax(111)
        phi = []
        tz0 = []
        tz0_th = []
        color = {10: 'tab:blue', 15: 'tab:orange', 20: 'tab:red'}
        for i, s in enumerate(self.stack):
            if type(s.data.data['phi_chamber']) == float:
                phi.append(s.data.data['phi_cons'])
            else:
                phi.append(np.float(s.data.data['phi_chamber'].replace(',', '.')))
            tz0.append(s.get_tz0()/s.data.fps)
            tz0_th.append(self.stefan[i].tz0)
            ax.plot(
                phi[-1], tz0[-1]/tz0_th[-1], 'o',
                color=color[s.data.data['Ta_cons']],
            )
        # print(phi, tz0)
        # p = np.polyfit(phi, tz0, 1)
        # f = np.poly1d(p)
        # phi_fit = np.linspace(np.min(phi), np.max(phi), 100)
        # ax.plot(phi_fit, f(phi_fit), '-.', color='tab:red')
        ax.xlabel(r'$\phi$ (%)')
        ax.ylabel(r'$t_{z_0}/t_{z_{0}, S_{te}}$')
        fig.show()

    def qa_phi(self, param='phi_int'):
        fig = pl.Figure()
        ax = fig.add_ax(111)
        phi = []
        qa = []
        qste = []
        for i, s in enumerate(self.stack):
            phi.append(np.float(s.data.data['phi_chamber'].replace(',', '.')))
            qa.append(s.get_qa(N=50))
            z0 = (s.geom.zb-s.geom.z0)/s.data.data['px_mm']
            qste.append(2.2*(0-s.data.Tnuc)/z0)

        ax.plot(phi, np.abs(qa), 'sk', mfc='None')
        p = np.polyfit(phi, np.abs(qa), 1)
        f = np.poly1d(p)
        phi_fit = np.linspace(np.min(phi), np.max(phi), 100)
        ax.plot(phi_fit, f(phi_fit), '-.', color='tab:red')
        ax.xlabel(r'$\phi$ (%)')
        ax.ylabel(r'$q_a (W/m^2)$')
        fig.show()

        fig = pl.Figure()
        ax = fig.add_ax(111)
        ax.plot(phi, np.abs(qa)/np.array(qste), 'sk', mfc='None')
        ax.xlabel(r'$\phi$ (%)')
        ax.ylabel(r'$q_a / q_{S_{te}}$')
        fig.show()

    def Ta_max_phi(self):
        """Plot Ta_max vs phi."""
        fig = pl.Figure()
        ax = fig.add_ax(111)
        phi = []
        Tamax = []
        for i, s in enumerate(self.stack):
            phi.append(np.float(s.data.data['phi_chamber'].replace(',', '.')))
            T = s.get_T_for_airflux()
            Tamax.append(np.max(T))

        ax.plot(phi, np.abs(Tamax), 'sk', mfc='None')
        ax.xlabel(r'$\phi$ (%)')
        ax.ylabel(r'$T_{a,max} (W/m^2)$')
        fig.show()


if __name__ == '__main__':
    date = [
        '15-10-2018', '15-10-2018', '15-10-2018', '15-10-2018',
        '17-10-2018', '17-10-2018', '17-10-2018', '17-10-2018', '17-10-2018',
        '22-10-2018', '22-10-2018',
        # '20-11-2018', '20-11-2018',	'20-11-2018',
        # '22-11-2018', '22-11-2018'
    ]
    serie = [
        1, 3, 5, 11,
        1, 3, 5, 7, 9,
        3, 7,
        # 5, 7, 9,
        # 3, 5
    ]
    # date = ['20-11-2017']*6
    # serie = range(1, 7)
    # date\= [5, 2, 5, 2]
    gest = Load()
    for i in range(len(date)):
        print(date[i], serie[i])
        gest.add_path(date[i], serie[i])

    # gest.kinetic()
    # gest.kinetic_adim_zO_tz0()
    # gest.tz0_repr()
    # gest.qa_phi()
    # gest.Ta_max_phi()
    gest.tip_angle()
