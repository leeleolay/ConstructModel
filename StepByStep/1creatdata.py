AtomID = {}

inputfilename = 'capsule_id.txt'
outputfilename = 'capsule_step1.data'

with open(inputfilename,'r') as f:
    ncount = 0
    for line in f.readlines():
        ncount = ncount + 1
        line = line.strip().split()
        AtomID[str(ncount)] = [line[0], line[1], line[2]]

xlo=ylo=zlo=-40
xhi=yhi=zhi=40

with open(outputfilename,"w") as f:
    f.write(
        f'''LAMMPS data file

{str(len(AtomID))} atoms
1 atom types 
0 bonds
1 bond types

{str(xlo)} {str(xhi)} xlo xhi
{str(ylo)} {str(yhi)} ylo yhi
{str(zlo)} {str(zhi)} zlo zhi

Masses

1 1

Atoms #full

'''
    )
    for key in AtomID:
        f.write(
            key + " 1 1 0.0 " + str(AtomID[key][0]) + " " + str(AtomID[key][1]) + " " + str(AtomID[key][2]) + "\n"
        )