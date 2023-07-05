import scipy.spatial as ss
import numpy as np
import re
AtomID = {}
BondID = {}
AngleID = {}
DihedralID = {}
Anglearea = {}
coords = []
Bondlenth = {}

inputfilename = 'capsule_step3.data'
outputfilename = 'capsule.data'

with open(inputfilename,'r') as f:
    ncount = 0
    for line in f.readlines():
        ncount = ncount + 1
        if re.search('.*xlo xhi',line):
            line = line.strip().split()
            xlo = float(line[0])
            xhi = float(line[1])
        elif re.search('.*ylo yhi',line):
            line = line.strip().split()
            ylo = float(line[0])
            yhi = float(line[1])
        elif re.search('.*zlo zhi',line):
            line = line.strip().split()
            zlo = float(line[0])
            zhi = float(line[1])
        elif re.search('.*atoms',line):
            line = line.strip().split()
            natom = int(line[0])
        elif re.search('.*bonds',line):
            line = line.strip().split()
            nbond = int(line[0])
        elif re.search('.*angles',line):
            line = line.strip().split()
            nangle = int(line[0])
        elif re.search('.*dihedrals',line):
            line = line.strip().split()
            ndihedral = int(line[0])
        elif re.search('Atoms',line):
            natomstart = ncount + 2
            natomend = natomstart + natom - 1
        elif re.search('Bonds',line):
            nbondstart = ncount + 2
            nbondend = nbondstart + nbond - 1
        elif re.search('Angles',line):
            nanglestart = ncount + 2
            nangleend = nanglestart + nangle - 1
        elif re.search('Dihedrals',line):
            ndihedralstart = ncount + 2
            ndihedralend = ndihedralstart + ndihedral - 1
with open(inputfilename,'r') as f:
    ncount = 0
    for line in f.readlines():
        ncount = ncount + 1
        if ncount >= natomstart and ncount <= natomend:
            line = line.strip().split()
            AtomID[line[0]] = [line[3], line[4], line[5]]
        if ncount >= nbondstart and ncount <= nbondend:
            line = line.strip().split()
            BondID[line[0]] = [line[2],line[3]]
            Bondlenth[line[0]] = [line[4]]
        if ncount >= nanglestart and ncount <= nangleend:
            line = line.strip().split()
            AngleID[line[0]] = [line[2],line[3],line[4]]
            Anglearea[line[0]] = [line[5]]
        if ncount >= ndihedralstart and ncount <= ndihedralend:
            line = line.strip().split()
            DihedralID[line[0]] = [line[2],line[3],line[4],line[5]]

with open(outputfilename, "w") as f:
    f.write(
        f'''LAMMPS data file

{str(len(AtomID))} atoms
1 atom types 
{str(len(BondID))} bonds
1 bond types
{str(len(AngleID))} angles
1 angle types 
{str(len(DihedralID))} dihedrals
1 dihedral types

{str(xlo)} {str(xhi)} xlo xhi
{str(ylo)} {str(yhi)} ylo yhi
{str(zlo)} {str(zhi)} zlo zhi

Masses

1 1

Atoms

'''
    )
    for key in AtomID:
        f.write(
            key + " 1001 1 " + str(AtomID[key][0]) + " " + str(AtomID[key][1]) + " " + str(-float(AtomID[key][2])) + "\n"
        )
    f.write(
        '''
Bonds

'''
    )
    for key in BondID:
        f.write(
            key + " 1 " + str(BondID[key][0]) + " " + str(BondID[key][1]) + " " + str(Bondlenth[key][0]) + "\n"
        )
    f.write(
        '''
Angles

'''
    )
    for key in AngleID:
        f.write(
            key + " 1 " + str(AngleID[key][0]) + " " + str(AngleID[key][1]) + " " + str(AngleID[key][2]) + " " + str(Anglearea[key][0]) + "\n"
        )
    f.write(
        '''
Dihedrals
    
'''
    )
    for key in DihedralID:
        f.write(
            key + " 1 " + str(DihedralID[key][2]) + " " + str(DihedralID[key][1]) + " " + str(DihedralID[key][0]) + " " + str(DihedralID[key][3]) + "\n"
        )

