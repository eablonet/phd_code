import numpy as np
from numpy import sqrt, argmin
import matplotlib.pyplot as plt


mu = 1.5e-5
L = 3e-3
z0 = L
V = 0.75
phi = .4

ki = 2.22
kl = 0.55
ka = 25e-3

Di = 1.18e-6
Da = 21e-6
Dl = 1.31e-7

cpi = 2050
cpa = 1021
cpl = 4200

rhoi = 916
rhoa = 1.2
rhol = 999

Ta_inf = 20
Tm = 0
Tw = -15

Lf = 3.35e5

Ste = cpi*(Tm - Tw)/Lf
Stea = cpa*(Ta_inf - Tm)/Lf
Stel = cpl*(Ta_inf - Tm)/Lf
delta = np.linspace(4.7e-5, 1e-3, 2)

print('delta : ', delta)

dt = 1e-2
tz0_th = []
tz0 = []
for d in delta:
    z_appro = [0]
    z = [0]
    t = [0]
    print(delta)
    for i in range(int(1e5)):
        t.append(t[-1]+dt)
        # zi = sqrt(
        #     z_appro[i]**2 + 2*dt*(
        #         Di*Ste -
        #         rhoa/rhoi * Da/delta * Stea * z_appro[i]
        #     )
        # )
        # z_appro.append(zi)
        zi = sqrt(
                z[i]**2 + 2*dt*(
                    Di*Ste -
                    rhol/rhoi * Dl * Stel * z[i] / (
                        z0 - z[i] + kl/ka*d
                    )
                )
            )
        z.append(zi)

    # tz0_appro t[argmin(np.abs(z0 - np.array(z_appro)))]
    tz0.append(t[argmin(np.abs(z0 - np.array(z)))])
    tz0_th.append(t[argmin(abs(z0 - sqrt(2*Di*Ste)*sqrt(t)))])

plt.figure()
plt.plot(ka*(Ta_inf - Tm)/delta, np.array(tz0)/tz0_th[0], '--k')
plt.xlabel('qa(z0)')
plt.ylabel('tz0/tz0_th')
plt.show()


phi = [
    41.7, 42.7, 45, 54.5, 55.7,
    # 72
]
tz0 = [1.25, 1.36, 1.33, 1.68, 1.88, 3.3]
q = [
    4100, 5300, 4900, 8100, 9030,
    # 0
]

plt.figure()
plt.plot(phi, q, '.k')
plt.xlabel(r'$\phi$')
plt.ylabel(r'$q_a$')
plt.grid(True)


z = np.polyfit(phi, q, 1)
p = np.poly1d(z)
print(z)
plt.plot(phi, p(phi), '--r')
plt.show()
"""

plt.figure(figsize=(8, 4.5))
plt.plot(t, z_appro, ls='--', color='tab:orange')
plt.plot(t, z, ls='--', color='tab:red')
plt.plot(t, sqrt(2*Di*Ste)*sqrt(t), ls='--', color='tab:blue')
plt.plot(
    [t[0], t[-1]],
    [L, L], '--k'
)
plt.plot(
    [tz0, tz0],
    [0, L], '--', color='tab:red'
)
plt.plot(
    [tz0_appro, tz0_appro],
    [0, L], '--', color='tab:orange'
)
plt.plot(
    [tz0_th, tz0_th],
    [0, L], '--', color='tab:blue'
)
plt.show()
"""
