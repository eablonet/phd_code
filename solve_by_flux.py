import numpy as np
from numpy import sqrt, argmin
import matplotlib.pyplot as plt


mu = 1.5e-5
z0 = 2e-3
V = 0.75
phi = .4

ki = 2.22
kl = 0.55
ka = 25e-3
ks = 401

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

qa = np.linspace(-15000, -2000, 200)
e = np.linspace(0, 2e-2, 50)

dt = .01
tz0 = []
tz0_sub = []
plt.figure()
tmax = 80
N = int(tmax/dt)

for q in range(len(qa)):
    zi_qa = [0]
    zi_sub = [0]
    t = [0]
    for i in range(N):
        t.append(t[i]+dt)
        zi_qa.append(
            np.sqrt(
                zi_qa[i]**2 +
                2*dt*(
                    Di*Ste +
                    qa[q]*zi_qa[i]/(rhoi*Lf)
                )
            )
        )

    if q % 10 == 0:
        plt.plot(
            t, zi_qa,
            ls='--', marker='^',
            markevery=100, ms=8, lw=2,
            color='tab:blue'
        )
    if np.max(zi_qa) < z0:
        print('==== Warning ====')
        print('Not enought point to solve')
    tz0.append(t[argmin(np.abs(z0 - np.array(zi_qa)))])

for delta in range(len(e)):
    zi_sub = [0]
    t = [0]
    for i in range(N):
        t.append(t[i]+dt)
        c = (
            -zi_sub[i]**2 -
            2*ki/ks*e[delta]*zi_sub[i] -
            2*Di*Ste*dt
        )
        Delta = 4*(ki/ks)**2*e[delta]**2 - 4*c
        zi_sub.append(
            -ki/ks*e[delta] + 1/2*np.sqrt(Delta)
        )

    if delta % 10 == 0:
        plt.plot(
            t, zi_sub,
            ls='--', marker='^',
            markevery=100, ms=8, lw=2,
            color='tab:red'
        )
    if np.max(zi_sub) < z0:
        print('==== Warning ====')
        print('Not enought point to solve')
    tz0_sub.append(t[argmin(np.abs(z0 - np.array(zi_sub)))])

plt.plot([0, 20], [z0, z0], color='tab:red')
tz0_th = t[argmin(abs(z0 - sqrt(2*Di*Ste)*sqrt(t)))]

plt.figure()
plt.plot(
    qa, np.array(tz0)/tz0_th,
    ls='--', marker='.', markevery=5,
    lw=2, ms=8,
    color='tab:blue'
)
qs = [-ki*(Tm-Tw)/(z0+ki*i/ks) for i in e]
plt.plot(
    qs, np.array(tz0_sub)/tz0_th,
    ls='--', marker='.', markevery=5,
    lw=2, ms=8,
    color='tab:red'
)
plt.xlabel(r'$q$ ($W/m^2$)')
plt.ylabel(r'$t_{z0}/t_{z0,th}$')

phi = [
    41.7, 42.7, 45, 54.5, 55.7, 72
]
tz0_ex = [1.25, 1.36, 1.33, 1.68, 1.88, 3.3]
q = []
for t in tz0_ex:
    q.append(
        qa[argmin(abs(np.array(tz0)/tz0_th - t))]
    )

# plot q vs phi
# -----------
plt.figure()
plt.plot(
    phi, q,
    ls='none', marker='.', markevery=1,
    lw=2, ms=8,
)
plt.xlabel(r'$\phi$ (%)')
plt.ylabel(r'$q_a$ ($W/m^2$)')
plt.grid(True)

z = np.polyfit(phi, q, 1)
p = np.poly1d(z)
print(z)
plt.plot(phi, p(phi), '--', color='tab:red')

# plot q/q_ste vs phi
# -----------
qste = ki*(Tm-Tw)/z0
plt.figure()
plt.plot(
    phi, [-i/qste for i in q],
    ls='none', marker='.', markevery=1,
    lw=2, ms=8,
)
plt.xlabel(r'$\phi$ (%)')
plt.ylabel(r'$q_a/q_{ste}$')
plt.grid(True)

z = np.polyfit(phi, [-i/qste for i in q], 1)
p = np.poly1d(z)
print(z)
plt.plot(phi, p(phi), '--', color='tab:red')

plt.show()
