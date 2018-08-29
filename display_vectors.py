import numpy as np
import matplotlib.pyplot as plt

size = int(input("Enter the number of points: "))
with open('HeI_0.dat') as f:
    array = [[float(x) for x in line.split()] for line in f]

z1 = []
h1 = []

for i in range(size):
    z1.append(array[i][0])
    h1.append(array[i][1])


with open('HeI_2500000.dat') as f:
    array = [[float(x) for x in line.split()] for line in f]

z2 = []
h2 = []

for i in range(size):
    z2.append(array[i][0])
    h2.append(array[i][1])

with open('H_250000000000.dat') as f:
    array = [[float(x) for x in line.split()] for line in f]

z3 = []
h3 = []

for i in range(size):
    z3.append(array[i][0])
    h3.append(array[i][1])

plt.plot(h1,z1, 'r--', h2, z2, 'g--', h3, z3, 'b--')
plt.show()
