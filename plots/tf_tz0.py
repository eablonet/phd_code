import numpy as np

from matplotlib import pyplot as plt

tz0 = np.array([18.4, 27.7, 18.9, 28.2, 22.7, 34.9, 43.1, 23.3, 35.5, 43.8])
tf = np.array([21.9, 32.9, 22.6, 33.6, 27.1, 41.6, 51.4, 27.7, 42.3, 52.2])

dt = tf - tz0
print(np.mean(dt), np.std(dt))
print(np.mean(dt)/tz0)
print(np.mean(dt/tz0))

plt.plot(tz0, tf, '--ok')
plt.show()
