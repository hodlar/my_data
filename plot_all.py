import numpy as np
import matplotlib.pyplot as plt
import math

size = int(input("Enter the number of points: "))
file0 = 'H_0.dat'
file1 = 'H_10000.dat'
file2 = 'H_25000.dat'
file3 = 'H_40000.dat'
field1 = '150G'
field2 = '380G'
field3 = '600G'

with open(file0) as f:
    array = [[float(x) for x in line.split()] for line in f]

z1 = []
h1 = []

for i in range(size):
    z1.append(math.log10(array[i][0]+1))
    h1.append(math.log10(array[i][1]))


with open(file1) as f:
    array = [[float(x) for x in line.split()] for line in f]

z2 = []
h2 = []

for i in range(size):
    z2.append(math.log10(array[i][0]))
    h2.append(math.log10(array[i][1]))


with open(file2) as f:
    array = [[float(x) for x in line.split()] for line in f]

z3 = []
h3 = []

for i in range(size):
    z3.append(math.log10(array[i][0]))
    h3.append(math.log10(array[i][1]))

with open(file3) as f:
    array = [[float(x) for x in line.split()] for line in f]

z4 = []
h4 = []

for i in range(size):
    z4.append(math.log10(array[i][0]))
    h4.append(math.log10(array[i][1]))

plt.figure(num='densidades')
plt.xlabel('T')
plt.ylabel('Density')
plt.plot(z1,h1, 'm--', label='0G')
plt.plot(z2, h2, 'r--', label = field1)
plt.plot(z3, h3, 'g--', label = field2)
plt.plot(z4, h4, 'b--', label = field3)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0)



with open(file0) as f:
    array = [[float(x) for x in line.split()] for line in f]

z1 = []
h1 = []
zeros = []

for i in range(size):
    z1.append(array[i][0])
    h1.append(array[i][1])
    zeros.append(0)


with open(file1) as f:
    array = [[float(x) for x in line.split()] for line in f]

z2 = []
h2 = []

for i in range(size):
    z2.append(array[i][0])
    h2.append(array[i][1])


with open(file2) as f:
    array = [[float(x) for x in line.split()] for line in f]

z3 = []
h3 = []

for i in range(size):
    z3.append(array[i][0])
    h3.append(array[i][1])

with open(file3) as f:
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
plt.plot(z,rel1, 'm--', label='0 vs '+field1)
plt.plot(z, rel2, 'r--', label=field1 +'vs'+field2)
plt.plot(z, rel3, 'g--', label=field2 +'vs' +field3)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0)


with open(file0) as f:
    array = [[float(x) for x in line.split()] for line in f]

z1 = []
h1 = []
zeros = []

for i in range(size):
    z1.append(array[i][0])
    h1.append(array[i][1])
    zeros.append(0)


with open(file1) as f:
    array = [[float(x) for x in line.split()] for line in f]

z2 = []
h2 = []

for i in range(size):
    z2.append(array[i][0])
    h2.append(array[i][1])


with open(file2) as f:
    array = [[float(x) for x in line.split()] for line in f]

z3 = []
h3 = []

for i in range(size):
    z3.append(array[i][0])
    h3.append(array[i][1])

with open(file3) as f:
    array = [[float(x) for x in line.split()] for line in f]

z4 = []
h4 = []

for i in range(size):
    z4.append(array[i][0])
    h4.append(array[i][1])

zeros = np.subtract(h1,h1)
z = np.subtract(z1,zeros)
rel1 = np.subtract(h2,h1)
rel2 = np.subtract(h3,h1)
rel3 = np.subtract(h4,h1)

#rel2 = list(set(h3) - set(h2))
#rel3 = list(set(h4) - set(h3))

plt.figure(num='absoluto')
plt.xlabel('T')
plt.ylabel('Density')
plt.plot(z, rel1, 'm--', label=field1)
plt.plot(z, rel2, 'r--', label=field2)
plt.plot(z, rel3, 'g--', label=field3)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0)

plt.show()
