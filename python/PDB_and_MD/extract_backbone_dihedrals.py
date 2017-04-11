# usage: python extract_backbone_dihedrals.py -i input.psf -o dihedral_numbers.txt

"""

OVERVIEW:

This script extracts the index numbers of the backbone dihedrals in a particular
PSF file from CHARMM.

"""

# Written by Thomas Gurry (thomasgurry@gmail.com)
# Updated: 10/16/2012

import os, string
import numpy as np
from optparse import OptionParser

# Read in arguments for the script

usage = "%prog -i INPUT_PSF -o OUTPUT_TEXTFILE"
parser = OptionParser(usage)
parser.add_option("-i", "--input_psf", type="string", dest="inputpsf")
parser.add_option("-o", "--output_dihedrals", type="string", dest="output_dihedrals")
(options, args) = parser.parse_args()

if( not options.inputpsf ):
    parser.error("No PSF file specified.")
if( not options.output_dihedrals ):
    parser.error("No dihedral output file specified.")

# File I/O (open)

psffile = open(options.inputpsf,'r')
outputfile = open(options.output_dihedrals,'w')

# Find relevant line indices

all_lines = psffile.readlines()

for i in range(len(all_lines)):
    line = all_lines[i].split()
    if (len(line) > 1):
        if (line[1] == "!NATOM"):
            beginning_ind = i+1
        elif (line[1] == "!NBOND:"):
            end_ind = i-2
        elif (line[1] == "!NPHI:"):
            dihedral_beginning_ind = i+1
        elif (line[1] == "!NIMPHI:"):
            dihedral_end_ind = i-2

# Read in relevant dihedral atom numbers

line = all_lines[end_ind].split()
nresidues = int(line[2])
resnums = range(2,nresidues)
phi_nums = np.zeros(shape=(4,len(resnums)), dtype=np.int)
psi_nums = np.zeros(shape=(4,len(resnums)), dtype=np.int)
phi_angle_nums = [0]*len(resnums)
psi_angle_nums = [0]*len(resnums)
counter = 0

for j in resnums:
    for i in range(beginning_ind,end_ind+1):
        line = all_lines[i].split()
        # Phi angle
        if (line[2] == str(j-1) and line[4] == "C"):
            phi_nums[0,counter] = int(line[0])
        if (line[2] == str(j) and line[4] == "N"):
            phi_nums[1,counter] = int(line[0])
        if (line[2] == str(j) and line[4] == "CA"):
            phi_nums[2,counter] = int(line[0])
        if (line[2] == str(j) and line[4] == "C"):
            phi_nums[3,counter] = int(line[0])
        # Psi angle
        if (line[2] == str(j) and line[4] == "N"):
            psi_nums[0,counter] = int(line[0])
        if (line[2] == str(j) and line[4] == "CA"):
            psi_nums[1,counter] = int(line[0])
        if (line[2] == str(j) and line[4] == "C"):
            psi_nums[2,counter] = int(line[0])
        if (line[2] == str(j+1) and line[4] == "N"):
            psi_nums[3,counter] = int(line[0])
    counter = counter + 1

# Find dihedral angle number for each quadruplet

for j in range(len(resnums)):
    counter = 0
    for i in range(dihedral_beginning_ind, dihedral_end_ind+1):
        line = all_lines[i].split()
        if(line[0] == str(phi_nums[0,j]) and line[1] == str(phi_nums[1,j]) and line[2] == str(phi_nums[2,j]) and line[3] == str(phi_nums[3,j])):
            phi_angle_nums[j] = 2*(counter+1)-1
        if(line[4] == str(phi_nums[0,j]) and line[5] == str(phi_nums[1,j]) and line[6] == str(phi_nums[2,j]) and line[7] == str(phi_nums[3,j])):
            phi_angle_nums[j] = 2*(counter+1)
        if(line[0] == str(psi_nums[0,j]) and line[1] == str(psi_nums[1,j]) and line[2] == str(psi_nums[2,j]) and line[3] == str(psi_nums[3,j])):
            psi_angle_nums[j] = 2*(counter+1)-1
        if(line[4] == str(psi_nums[0,j]) and line[5] == str(psi_nums[1,j]) and line[6] == str(psi_nums[2,j]) and line[7] == str(psi_nums[3,j])):
            psi_angle_nums[j] = 2*(counter+1)
        counter = counter + 1



teststr1 = ''
teststr2 = ''

for i in range(len(phi_angle_nums)):
    teststr1 = teststr1+' '+str(phi_angle_nums[i])

for i in range(len(psi_angle_nums)):   
    teststr2 = teststr2+' '+str(psi_angle_nums[i])

outputfile.write(teststr1+'\n')
outputfile.write(teststr2+'\n')

# File I/O (close)

psffile.close()
outputfile.close()
