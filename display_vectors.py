import numpy as np
import matplotlib.pyplot as plt

size = int(input("Enter the number of points: "))
with open('hydrostatic.dat') as f:
    array = [[float(x) for x in line.split()] for line in f]

z1 = []
h1 = []

for i in range(size):
    z1.append(array[i][1])
    h1.append(array[i][4])


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

plt.figure(num='densidades')
plt.xlabel('T')
plt.ylabel('Density')
plt.plot(z1,h1, 'm--', z2, h2, 'r--', z3, h3, 'g--', z4, h4, 'b--')
#plt.plot(z1,h1, 'r--')
plt.show()
