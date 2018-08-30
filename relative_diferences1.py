import numpy as np
import matplotlib.pyplot as plt

size = int(input("Enter the number of points: "))
with open('hydrostatic.dat') as f:
    array = [[float(x) for x in line.split()] for line in f]

z1 = []
h1 = []
zeros = []

for i in range(size):
    z1.append(array[i][1])
    h1.append(array[i][4])
    zeros.append(0)


with open('H_25000.dat') as f:
    array = [[float(x) for x in line.split()] for line in f]

z2 = []
h2 = []

for i in range(size):
    z2.append(array[i][0])
    h2.append(array[i][1])


with open('H_2500000.dat') as f:
    array = [[float(x) for x in line.split()] for line in f]

z3 = []
h3 = []

for i in range(size):
    z3.append(array[i][0])
    h3.append(array[i][1])

with open('H_250000000000.dat') as f:
    array = [[float(x) for x in line.split()] for line in f]

z4 = []
h4 = []

for i in range(size):
    z4.append(array[i][0])
    h4.append(array[i][1])

zeros = np.subtract(h1,h1)
z = np.subtract(z1,zeros)
rel1 = np.subtract(h2,h1)
rel2 = np.subtract(h3,h2)
rel3 = np.subtract(h4,h3)


plt.figure(num='relativo')
plt.xlabel('T')
plt.ylabel('Density')
plt.plot(z,rel1, 'm--', z, rel2, 'r--', z, rel3, 'g--')
#plt.plot(z,rel1, 'r--')
plt.show()
