import scipy.interpolate as interp
import matplotlib.pyplot as plt
import numpy as np
import hermiteSpline as hs
from IPython import embed


plt.title("Centroidal Voronoi Tesselation Energy", size=16)
plt.xlabel("Number of Generators", size=12)
plt.ylabel("Energy (1e9)", size=12)
plt.xlim(1, 17)
#plt.ylim(-100000000,7000000000)

y = np.array((6463169658.5099, 1271863820.59449, 472534881.328941,
              263864995.336371,132889255.869661, 60904655.442872,
              35020018.020686, 24719313.79766, 16105092.21555,
              12296504.406197, 11132316.653403, 10136002.913738,
              9626162.703109, 9515755.897589, 9498464.90249))

y /= 1000000000

plt.ylim(-.1,7)

x = np.array((2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16))

sp = hs.HermiteSpline(np.array((x, y)))
spx = np.arange(0, sp.length, sp.length / 1000)
spy = sp(spx)

f = interp.interp1d(x, y, 5)
xNew = np.arange(2, 16, 0.1)
yNew = f(xNew)

#embed()
#plt.plot(spy)
plt.plot(xNew,yNew,linewidth=3)
plt.show()
