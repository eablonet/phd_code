import packages.gui.ploter as plt
tnuc = [
    80,
    223,
    126,
    147,
    283,
    286,
    0,
    14,
    0,
]
tf = [
    40,
    43,
    51,
    46,
    46,
    38,
    80,
    76,
    70,
]

fig = plt.Figure()
a = fig.add_ax(111)
a.plot(range(len(tnuc)), tnuc, marker='s', mfc='none', color='tab:blue', label='Temps de nucléation')
a.plot(range(len(tnuc)), tf,  marker='s', mfc='none', color='tab:red', label='Temps de solidification')
a.xlabel('Serie')
a.ylabel('Temps (s)')
a.legend()
fig.show()
