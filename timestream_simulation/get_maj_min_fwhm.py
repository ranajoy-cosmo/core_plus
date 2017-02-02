import numpy as np
import sys

mid = np.float(sys.argv[1])
elp = np.float(sys.argv[2])

x = ((2 + elp)/(2 - elp))**2

print x

maj = np.sqrt(2)*mid*x/np.sqrt(1 + x)
min = np.sqrt(2)*mid/np.sqrt(1 + x)

print maj, min
