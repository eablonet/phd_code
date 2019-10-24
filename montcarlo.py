# monte-carlo pi calculation

import numpy as np
d = 8
N = int(1e5)
R = []
for i in range(N):
    x = np.random.rand(d, 1)
    R.append(np.sum(np.power(x, 2)))

R = np.array(R)
p = 2**d*24 * len(R[R <= 1]) / N
print('pi = ', p**(1/4))
