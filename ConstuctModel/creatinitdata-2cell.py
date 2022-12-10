import scipy.spatial as ss
import numpy as np
import re

cAtomID = {}
cBondID = {}
cAngleID = {}
cDihedralID = {}
cAnglearea = {}
ccoords = []
cBondlenth = {}

bAtomID = {}
bBondID = {}
bAngleID = {}
bDihedralID = {}
bAnglearea = {}
bcoords = []
bBondlenth = {}

inputfilename1 = 'capsule.data'
inputfilename2 = 'ball.data'
outputfilename = 'init.data'

xlo = -25
xhi = 25
ylo = -12.5
yhi = 12.5
zlo = -150
zhi = 150

particle_xlo = -25
particle_xhi = 25
particle_ylo = -12.5
particle_yhi = 12.5
particle_zlo = 130
particle_zhi = 138
disDropAtom = 3
numberofparticle = np.arange(particle_xlo,particle_xhi,disDropAtom).size * np.arange(particle_ylo,particle_yhi,disDropAtom) .size * np.arange(particle_zlo,particle_zhi,disDropAtom).size
distBetweenCellAndBall = -50
distOfCell = 5 + 20.0
distOfball = 15.0
shiftx = -12.5

cellxnum = 2
cellynum = 1
numOfCell = cellxnum*cellynum
ballxnum = 0 #14
ballynum = 0 #14
ballznum = 0 #1
numOfBall = ballxnum*ballynum*ballznum

distwallupper = zhi - 10
distwalldown = zlo + 10
distbetweenparticleofwall = 0.5
deltaedge = 0.1
numparticlefromwall = len(np.arange(xlo+deltaedge,xhi-deltaedge,distbetweenparticleofwall))*len(np.arange(ylo+deltaedge,yhi-deltaedge,distbetweenparticleofwall)) 

with open(inputfilename,'r') as f:
    ncount = 0
    for line in f.readlines():
        ncount = ncount + 1
        if re.search('.*xlo xhi',line):
            line = line.strip().split()
            cxlo = float(line[0])
            cxhi = float(line[1])
        elif re.search('.*ylo yhi',line):
            line = line.strip().split()
            cylo = float(line[0])
            cyhi = float(line[1])
        elif re.search('.*zlo zhi',line):
            line = line.strip().split()
            czlo = float(line[0])
            czhi = float(line[1])
        elif re.search('.*atoms',line):
            line = line.strip().split()
            cnatom = int(line[0])
        elif re.search('.*bonds',line):
            line = line.strip().split()
            cnbond = int(line[0])
        elif re.search('.*angles',line):
            line = line.strip().split()
            cnangle = int(line[0])
        elif re.search('.*dihedrals',line):
            line = line.strip().split()
            cndihedral = int(line[0])
        elif re.search('Atoms',line):
            cnatomstart = ncount + 2
            cnatomend = cnatomstart + cnatom - 1
        elif re.search('Bonds',line):
            cnbondstart = ncount + 2
            cnbondend = cnbondstart + cnbond - 1
        elif re.search('Angles',line):
            cnanglestart = ncount + 2
            cnangleend = cnanglestart + cnangle - 1
        elif re.search('Dihedrals',line):
            cndihedralstart = ncount + 2
            cndihedralend = cndihedralstart + cndihedral - 1
with open(inputfilename,'r') as f:
    ncount = 0
    for line in f.readlines():
        ncount = ncount + 1
        if ncount >= cnatomstart and ncount <= cnatomend:
            line = line.strip().split()
            cAtomID[line[0]] = [line[3], line[4], line[5]]
        if ncount >= cnbondstart and ncount <= cnbondend:
            line = line.strip().split()
            cBondID[line[0]] = [line[2],line[3]]
            cBondlenth[line[0]] = [line[4]]
        if ncount >= cnanglestart and ncount <= cnangleend:
            line = line.strip().split()
            cAngleID[line[0]] = [line[2],line[3],line[4]]
            cAnglearea[line[0]] = [line[5]]
        if ncount >= cndihedralstart and ncount <= cndihedralend:
            line = line.strip().split()
            cDihedralID[line[0]] = [line[2],line[3],line[4],line[5]]

with open(inputfilename2,'r') as f:
    ncount = 0
    for line in f.readlines():
        ncount = ncount + 1
        if re.search('.*xlo xhi',line):
            line = line.strip().split()
            bxlo = float(line[0])
            bxhi = float(line[1])
        elif re.search('.*ylo yhi',line):
            line = line.strip().split()
            bylo = float(line[0])
            byhi = float(line[1])
        elif re.search('.*zlo zhi',line):
            line = line.strip().split()
            bzlo = float(line[0])
            bzhi = float(line[1])
        elif re.search('.*atoms',line):
            line = line.strip().split()
            bnatom = int(line[0])
        elif re.search('.*bonds',line):
            line = line.strip().split()
            bnbond = int(line[0])
        elif re.search('.*angles',line):
            line = line.strip().split()
            bnangle = int(line[0])
        elif re.search('.*dihedrals',line):
            line = line.strip().split()
            bndihedral = int(line[0])
        elif re.search('Atoms',line):
            bnatomstart = ncount + 2
            bnatomend = bnatomstart + bnatom - 1
        elif re.search('Bonds',line):
            bnbondstart = ncount + 2
            bnbondend = bnbondstart + bnbond - 1
        elif re.search('Angles',line):
            bnanglestart = ncount + 2
            bnangleend = bnanglestart + bnangle - 1
        elif re.search('Dihedrals',line):
            bndihedralstart = ncount + 2
            bndihedralend = bndihedralstart + bndihedral - 1
with open(inputfilename2,'r') as f:
    ncount = 0
    for line in f.readlines():
        ncount = ncount + 1
        if ncount >= bnatomstart and ncount <= bnatomend:
            line = line.strip().split()
            bAtomID[line[0]] = [line[3], line[4], line[5]]
        if ncount >= bnbondstart and ncount <= bnbondend:
            line = line.strip().split()
            bBondID[line[0]] = [line[2],line[3]]
            bBondlenth[line[0]] = [line[4]]
        if ncount >= bnanglestart and ncount <= bnangleend:
            line = line.strip().split()
            bAngleID[line[0]] = [line[2],line[3],line[4]]
            bAnglearea[line[0]] = [line[5]]
        if ncount >= bndihedralstart and ncount <= bndihedralend:
            line = line.strip().split()
            bDihedralID[line[0]] = [line[2],line[3],line[4],line[5]]

with open(outputfilename, "w") as f:
    f.write(
        f'''LAMMPS data file

{str(len(cAtomID)*numOfCell+len(bAtomID)*numOfBall+numparticlefromwall*2+numberofparticle)} atoms
3 atom types 
{str(len(cBondID)*numOfCell+len(bBondID)*numOfBall)} bonds
2 bond types
{str(len(cAngleID)*numOfCell+len(bAngleID)*numOfBall)} angles
2 angle types 
{str(len(cDihedralID)*numOfCell+len(bDihedralID)*numOfBall)} dihedrals
2 dihedral types

{str(xlo)} {str(xhi)} xlo xhi
{str(ylo)} {str(yhi)} ylo yhi
{str(zlo)} {str(zhi)} zlo zhi

Masses

1 1
2 1
3 1

Atoms

'''
    )
    # create cell
    for i in range(0,cellxnum):
        for j in range(0,cellynum):
            for key in cAtomID:
                f.write(
                    str(int(key)+len(cAtomID)*(i*cellynum+j)) + " " + str(1001+i*cellynum+j)+ " 1 " + str(float(cAtomID[key][0])+distOfCell*i+shiftx) + " " + str(float(cAtomID[key][1])+distOfCell*j) + " " + str(float(cAtomID[key][2])) + "\n"
                )
    # create droplet
    for i in range(0,ballxnum):
        for j in range(0,ballynum):
            for k in range(0,ballznum):
                for key in bAtomID:
                    f.write(
                        str(int(key)+len(cAtomID)*numOfCell+len(bAtomID)*(i*ballynum*ballznum+j*ballznum+k)) + " " + str(1001+numOfCell+(i*ballynum*ballznum+j*ballznum+k)) + " 2 " + str(float(bAtomID[key][0])+distOfball*i+0) + " " + str(float(bAtomID[key][1])+distOfball*j+0) + " " + str(float(bAtomID[key][2])-distOfball*k+distBetweenCellAndBall) + "\n"
                    )
    # create upper wall
    for i in np.arange(xlo+deltaedge,xhi-deltaedge,distbetweenparticleofwall):
        for j in np.arange(ylo+deltaedge,yhi-deltaedge,distbetweenparticleofwall):
            number = 1
            f.write(
                str(int(number)+len(cAtomID)*numOfCell+len(bAtomID)*numOfBall) + " " + str(1001+numOfCell+numOfBall+number) + " 2 " + str(i) + " " + str(j) + " " + str(distwallupper) + "\n"
            )
            number += 1    
    # create down wall
    for i in np.arange(xlo+deltaedge,xhi-deltaedge,distbetweenparticleofwall):
        for j in np.arange(ylo+deltaedge,yhi-deltaedge,distbetweenparticleofwall):
            number = 1
            f.write(
                str(int(number)+len(cAtomID)*numOfCell+len(bAtomID)*numOfBall+int(numparticlefromwall)) + " " + str(1001+numOfCell+numOfBall+number) + " 2 " + str(i) + " " + str(j) + " " + str(distwalldown) + "\n"
            )
            number += 1
    # create droplet atoms
    for i in np.arange(particle_xlo,particle_xhi,disDropAtom):
        for j in np.arange(particle_ylo,particle_yhi,disDropAtom):
            for k in np.arange(particle_zlo,particle_zhi,disDropAtom):
                number = 1
                f.write(
                    str(int(number)+len(cAtomID)*numOfCell+len(bAtomID)*numOfBall+int(numparticlefromwall)*2) + " " + str(1001+numOfCell+numOfBall+2+number) + " 3 " + str(i) + " " + str(j) + " " + str(k) + "\n"
                )
                number += 1
    f.write(
        '''
Bonds

'''
    )
    for i in range(numOfCell):
        for key in cBondID:
            f.write(
                str(int(key)+len(cBondID)*i) + " 1 " + str(int(cBondID[key][0])+len(cAtomID)*i) + " " + str(int(cBondID[key][1])+len(cAtomID)*i) + " " + str(cBondlenth[key][0]) + "\n"
            )
    for i in range(numOfBall):
        for key in bBondID:
            f.write(
                str(int(key)+len(bBondID)*i+len(cBondID)*numOfCell) + " 2 " + str(int(bBondID[key][0])+len(bAtomID)*i+len(cAtomID)*numOfCell) + " "  + str(int(bBondID[key][1])+len(bAtomID)*i+len(cAtomID)*numOfCell) + " " + str(bBondlenth[key][0]) + "\n"
            )
    f.write(
        '''
Angles

'''
    )
    for i in range(numOfCell):
        for key in cAngleID:
            f.write(
                str(int(key)+len(cAngleID)*i) + " 1 " + str(int(cAngleID[key][0])+len(cAtomID)*i) + " " + str(int(cAngleID[key][1])+len(cAtomID)*i) + " " + str(int(cAngleID[key][2])+len(cAtomID)*i) + " " + str(cAnglearea[key][0]) + "\n"
            )
    for i in range(numOfBall):
        for key in bAngleID:
            f.write(
                str(int(key)+len(bAngleID)*i+len(cAngleID)*numOfCell) + " 2 " + str(int(bAngleID[key][0])+len(bAtomID)*i+len(cAtomID)*numOfCell) + " " + str(int(bAngleID[key][1])+len(bAtomID)*i+len(cAtomID)*numOfCell) + " " + str(int(bAngleID[key][2])+len(bAtomID)*i+len(cAtomID)*numOfCell) + " " + str(bAnglearea[key][0]) + "\n"
            )
    f.write(
        '''
Dihedrals
    
'''
    )
    for i in range(numOfCell):
        for key in cDihedralID:
            f.write(
                str(int(key)+len(cDihedralID)*i) + " 1 " + str(int(cDihedralID[key][2])+len(cAtomID)*i) + " " + str(int(cDihedralID[key][1])+len(cAtomID)*i) + " " + str(int(cDihedralID[key][0])+len(cAtomID)*i) + " " + str(int(cDihedralID[key][3])+len(cAtomID)*i) + "\n"
            )
    for i in range(numOfBall):
        for key in bDihedralID:
            f.write(
                str(int(key)+len(bDihedralID)*i+len(cDihedralID)*numOfCell) + " 2 " + str(int(bDihedralID[key][2])+len(bAtomID)*i+len(cAtomID)*numOfCell) + " " + str(int(bDihedralID[key][1])+len(bAtomID)*i+len(cAtomID)*numOfCell) + " " + str(int(bDihedralID[key][0])+len(bAtomID)*i+len(cAtomID)*numOfCell) + " " + str(int(bDihedralID[key][3])+len(bAtomID)*i+len(cAtomID)*numOfCell) + "\n"
            )

