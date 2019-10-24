# Calculate the volume for a condensation flux.

import packages.gui.ploter as pl
import numpy as np

qa = np.array([15000, 10000, 7000, 5000, 2000])
r = 3e-3
h = 1e-3
t = 10
rho = 999
Lv = 2264e3
v = qa*np.pi*(r**2+h**2)*t / (rho*Lv)

fig = pl.Figure()
ax = fig.add_ax(111)
ax.plot(qa, v, color='tab:red', marker='s')
ax.xlabel('Condensation flux')
ax.ylabel('Volume')
fig.show()
