# usage: python calc_RG.py -i pdbdir -o output.dat

# Calculates Radius of Gyration from each PDB in a folder

import os, sys
from optparse import OptionParser
from math import sqrt

def CalcRG(infile):
    # Calculate C-alpha radius of gyration
    coords = []
    pdb = open(infile,"r")
    for line in pdb:
        l = line.split()
        if l[0] == 'ATOM':
            # add a space between coordinates to make sure line splits correctly
            L = line[0:38]+" "+line[38:46]+" "+line[46:]; 
            LL = L.split()
            if LL[2] == "CA":
                try:
                    float(LL[4])
                    coords.append( [float(x) for x in LL[5:8]] )
                except:
                    coords.append( [float(x) for x in LL[6:9]] )
    pdb.close()
    n = len(coords)
    cm = [sum([x[i] for x in coords])/n for i in range(3)]
    R = sum([sum([coords[i][j]**2 for j in range(3)]) for i in range(len(coords))])/n
    r = sum([cm[i]**2 for i in range(3)])
    return sqrt(R-r)

# Read in arguments
usage = "%prog -i pdbdir -o output.dat"
parser = OptionParser(usage)
parser.add_option("-i", "--pdbdir", dest="pdbdir", help="path to directory containing PDBs")
parser.add_option("-o", "--outfile", dest="outfile", help="output file")
(options, args) = parser.parse_args()

if( not options.pdbdir ):
    parser.error("No PDB directory specified.")
if (not options.outfile ):
    parser.error("No output file specified.")

pdbdir = options.pdbdir
outfile = options.outfile

pdbfiles = os.listdir(pdbdir)
rgvals = len(pdbfiles)*[0]

for i in range(len(pdbfiles)):
    rgvals[i] = CalcRG(pdbdir+pdbfiles[i])

meanRG = sum(rgvals)/len(rgvals)

fout = open(outfile,'w')
fout.write("Ensemble average radius of gyration: " + str(meanRG) + " Angstroms")
fout.close()
