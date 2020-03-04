'''
Plots variables from log files in current directory
'''

from matplotlib import pyplot as plt
import numpy as np
import os

path = os.getcwd()
data_names = []

for root, dirs, files in os.walk(path, topdown=True):
    for name in files:
        if name[-18:] == 'thermalization.log':
            #val = [a for i in file.split('_') for a in i.split('-') if a[0]=='tstep']
            val = (root.split('/')[-1], name.split('.log')
                   [0], os.path.join(root, name))
            data_names.append(val)

print(len(data_names))

titles = set([i[0] for i in data_names])
subtitle = set([i[1].split('_')[-5] for i in data_names])
print(titles)
print(subtitle)
for t in titles:
    for st in subtitle:
        fig1 = plt.figure(1)
        ax1_1 = fig1.add_subplot(1, 1, 1)
        ax1_1.set_title('{}\n{}'.format(t, st))
        for d in data_names:
            if d[0] == t and d[1].split('_')[-5] == st:
                label = d[1].split('_')[-3]+'_'+d[1].split('_')[-2]
                data = np.genfromtxt(fname=d[2], skip_header=True)
                ax1_1.plot(data[:, 0], data[:, 8], label=label)

        ax1_1.set_ylabel('Temperature / -')
        ax1_1.set_xlabel('Time Step / -')
        ax1_1.legend()
        fig1.tight_layout()
        fig1.savefig('./Temperature_{}_{}.png'.format(t, st))
        plt.clf()


# ax1_1.xaxis.set_major_locator(plt.MaxNLocator(8))
