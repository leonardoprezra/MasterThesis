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

    new_t.append(t[frames])

# Convert .gsd to .pos


for in_file in args.in_file:
    # Output .gsd file name
    outFN = dir_name + \
        in_file.split('/')[-1][:-4]+'_frame-{}.gsd'.format(args.frame_num)

    slicegsd(frames=args.frame_num, inFN=in_file, outFN=outFN)
