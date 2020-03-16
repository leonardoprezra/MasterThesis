'''
Extract the last <frame_num> of the <in_file> and saves them as .gsd and .pos.

Positional arguments:

frame_num : number of frames to extract from the end of the .gsd file
in_file : input .gsd file
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
parser.add_argument('frame_num', type=int,
                    help='number of frames to extract from the end of the .gsd file')
parser.add_argument('in_file', help='input .gsd file')

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

    new_t.extend(t[-10:])

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

    with open(inFN[:-3] + 'pos', 'w') as file1:

        for f in t:

            box = f.configuration.box
            file1.write("box {:f} {:f} {:f}\n".format(box[0], box[1], box[2]))

            for i in range(f.particles.N):
                if f.particles.typeid[i] == 0:
                    file1.write("sphere {:f} {:s} ".format(f.particles.diameter[i], 'FF0000') + "{:f} {:f} {:f}\n".format(
                        *f.particles.position[i]))

                elif f.particles.typeid[i] == 1:
                    file1.write("sphere {:f} {:s} ".format(f.particles.diameter[i], '0000FF') + "{:f} {:f} {:f}\n".format(
                        *f.particles.position[i]))

            file1.write("eof\n")


# Output .gsd file name
outFN = dir_name + args.in_file.split('/')[-1]

slicegsd(frames=args.frame_num, inFN=args.in_file, outFN=outFN)

convert2pos(inFN=outFN)
