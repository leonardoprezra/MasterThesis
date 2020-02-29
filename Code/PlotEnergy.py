from matplotlib import pyplot as plt
import numpy as np
import os

path = '/nishome/students/leonardo/Dokumente/Thesis/Code/try'
data_names = []

for root, dirs, files in os.walk(path, topdown=True):
    for name in files:
        if name[-3:] == 'log':
            #val = [a for i in file.split('_') for a in i.split('-') if a[0]=='tstep']
            val = (root.split('/')[-1], name.split('.log')
                   [0], os.path.join(root, name))
            data_names.append(val)

print(len(data_names))

titles = set([i[0] for i in data_names])

for t in titles:
    fig1 = plt.figure(1)
    ax1_1 = fig1.add_subplot(1, 1, 1)
    ax1_1.set_title(t)
    for d in data_names:
        if d[0] == t:
            label = d[1]
            data = np.genfromtxt(fname=d[2], skip_header=True)
            ax1_1.plot(data[:, 0], data[:, 4], label=label)

    ax1_1.set_ylabel('Potential Energy / -')
    ax1_1.set_xlabel('Time / -')
    ax1_1.legend()
    fig1.tight_layout()
    fig1.savefig('./PotEnergy_{}.png'.format(t))
    plt.clf()


# ax1_1.xaxis.set_major_locator(plt.MaxNLocator(8))
