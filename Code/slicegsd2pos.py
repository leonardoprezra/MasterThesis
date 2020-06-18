'''
Extract the last <frame_num> of the <in_file> and saves them as .gsd and .pos.

'''

import gsd
import gsd.hoomd
import math

import sys
import argparse
import os

# Name of directory to store converted files
dir_name = 'data_slice_pos/'

# Command line argument parsing
parser = argparse.ArgumentParser(
    description='Extract the last <frame_num> of the <in_file> and saves them as .gsd and .pos.')
parser.add_argument('-n', '--frame-number', type=int, dest='frame_num',
                    help='number of frames to extract from the end of the .gsd file')
parser.add_argument('-i', '--in-file', type=str, nargs='+',
                    dest='in_file', help='input .gsd file')

args = parser.parse_args()

# Create directories
try:
    if(not os.path.exists(dir_name)):
        os.mkdir(dir_name)
except OSError as e:
    if e.errno != 17:
        raise
    pass


# Extract frames from .gsd file


def slicegsd(frames, inFN, outFN):
    '''Extract the last <frames> from <inFN> and writes them in a new file.

    Parameters
    ----------
    frames : int
        Number of frames to extract from the end of a .gsd file.
    inFN : str
        Input .gsd file.
    outFN : str
        Output .gsd file.

    Returns
    -------
    None
    '''
    t = gsd.hoomd.open(name=inFN, mode='rb')

    new_t = gsd.hoomd.open(name=outFN, mode='wb')

    new_t.extend(t[-frames:])

# Convert .gsd to .pos


def convert2pos(inFN):
    '''Convert .gds file to .pos.

    Parameters
    ----------
    inFN : str
        Input .gsd file.

    Returns
    -------
    None
    '''
    t = gsd.hoomd.open(name=inFN, mode='rb')
    print('!!!!!!!!!!!!!!!!!!!')
    print(inFN)
    print(inFN.split('/')[1].split('_')[1].split('-')[1])

    shape = inFN.split('/')[1].split('_')[1].split('-')[1]
    dimensions = int(inFN.split('_')[4].split('-')[1])

    if shape != 'one':
        # Diameter of core particles
        for i in range(t[0].particles.N):
            if t[0].particles.typeid[i] == 0:
                diam_core = t[0].particles.diameter[i]
                break

        # Diameter of halo particles
        for i in range(t[0].particles.N):
            if t[0].particles.typeid[i] == 1:
                diam_halo = t[0].particles.diameter[i]
                break
    else:
        # Diameter of core particles
        for i in range(t[0].particles.N):
            if t[0].particles.typeid[i] == 0:
                diam_core = t[0].particles.diameter[i]
                break

    with open(inFN[:-3] + 'pos', 'w') as file1:
        if dimensions == 2:
            for f in t:
                box = f.configuration.box

                if shape != 'one':
                    file1.write("box {:f} {:f} {:f}\n".format(
                        box[0], box[1], 0))
                    file1.write(
                        'def A "sphere {:f} FF0000"\n'.format(diam_core))  # Red
                    file1.write('def B "sphere {:f} 0000FF"\n'.format(
                        diam_halo))  # Blue
                else:
                    file1.write("box {:f} {:f} {:f}\n".format(
                        box[0], box[1], 0))
                    file1.write(
                        'def A "sphere {:f} FF0000"\n'.format(diam_core))  # Red

                for i in range(f.particles.N):
                    if f.particles.typeid[i] == 0:
                        file1.write("A " + "{:f} {:f} {:f}\n".format(
                            *f.particles.position[i]))

                    elif f.particles.typeid[i] == 1:
                        file1.write("B " + "{:f} {:f} {:f}\n".format(
                            *f.particles.position[i]))

                file1.write("eof\n")

        elif dimensions == 3:
            for f in t:

                box = f.configuration.box

                if shape != 'one':
                    file1.write("box {:f} {:f} {:f}\n".format(
                        box[0], box[1], box[2]))
                    file1.write(
                        'def A "sphere {:f} FF0000"\n'.format(diam_core))  # Red
                    file1.write('def B "sphere {:f} 0000FF"\n'.format(
                        diam_halo))  # Blue
                else:
                    file1.write("box {:f} {:f} {:f}\n".format(
                        box[0], box[1], box[2]))
                    file1.write(
                        'def A "sphere {:f} FF0000"\n'.format(diam_core))  # Red

                for i in range(f.particles.N):
                    if f.particles.typeid[i] == 0:
                        file1.write("A " + "{:f} {:f} {:f}\n".format(
                            *f.particles.position[i]))

                    elif f.particles.typeid[i] == 1:
                        file1.write("B " + "{:f} {:f} {:f}\n".format(
                            *f.particles.position[i]))

                file1.write("eof\n")


for in_file in args.in_file:
    # Output .gsd file name
    outFN = dir_name + in_file.split('/')[-1]

    slicegsd(frames=args.frame_num, inFN=in_file, outFN=outFN)

    convert2pos(inFN=outFN)
