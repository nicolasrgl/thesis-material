import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def calculate_average(x_axis, y_axis):
    avg = 0.0

    for i in range(0,len(x_axis)):
        avg += x_axis[i] * y_axis[i]
    
    return (avg / 1000.0)


file = open(sys.argv[1], 'r')
x_axis = []
y_axis = []
x_axis2 = []
y_axis2 = []

while 1:
    line = file.readline()
    if(line.startswith('-')):
        break
    if(line.startswith('0')):
        continue
    line = line.split()
    x_axis.append(float(line[0])-0.2)
    y_axis.append(int(line[1]))

fig, ax = plt.subplots()
ax.bar(x_axis, y_axis, width=0.4, color='red')

while 1:
    line = file.readline()
    if(line.startswith('-')):
        break
    if(line.startswith('0')):
        continue
    line = line.split()
    x_axis2.append(float(line[0])+0.2)
    y_axis2.append(int(line[1]))

dataset = sys.argv[1].split('/')[1].split('.')[0].split('_')
    
ax.bar(x_axis2, y_axis2, width=0.4, color='blue')

ax.set_ylabel('Results for given size')
ax.set_xticks(np.arange(min(x_axis2)-0.2, max(x_axis2)+0.8, max(x_axis2)//25))
ax.set_xlabel('Size of Certificate')
ax.set_title(dataset[1].upper() + ' certificate size comparison on ' + dataset[0] + ' curve set')
plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
fig.tight_layout()

sc = round(calculate_average(x_axis, y_axis), 4)
fl = round(calculate_average(x_axis2, y_axis2), 4)

print("Average certificate length SHORTEST_CERTIFICATE: " + str(sc))
print("Average certificate length FRECHET_LIGHT: " + str(fl) + "\n")

print("Improvement factor: " + str(fl / sc))

plt.show()
