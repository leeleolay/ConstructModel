from collections import Counter
import networkx as nx
import re
import math
import scipy.spatial as ss
import numpy as np

AtomID = {}
Atombond = {}
BondID = {}
Threeatom = []
AngleID = {}
DihedralID = {}
bondlenth = {} 
anglearea = {}

inputfilename1 = 'capsule_step2.data'
inputfilename2 = 'capsule_angle.txt'
outputfilename = 'capsule_step3.data'

class Point:
    def __init__(self,x=0,y=0,z=0):
        self.x=x
        self.y=y
        self.z=z
    def getx(self):
        return self.x
    def gety(self):
        return self.y
    def getz(self):
        return self.z
class Getlen:
    def __init__(self,p1,p2):
        self.x=p1.getx()-p2.getx()
        self.y=p1.gety()-p2.gety()
        self.z=p1.getz()-p2.getz()
        self.len = math.sqrt((self.x**2)+(self.y**2)+(self.z**2))
    def getlen(self):
        return self.len
class Getarea:
    def __init__(self,l1,l2,l3):
        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.p = (self.l1+self.l2+self.l3)/2
        self.area = math.sqrt(self.p*(self.p-self.l1)*(self.p-self.l2)*(self.p-self.l3))
    def getarea(self):
        return self.area

def insert_sort(list): # sort of order list
    n = len(list)
    for i in range(1, n):
        # 后一个元素和前一个元素比较
        # 如果比前一个小
        if list[i] < list[i - 1]:
            # 将这个数取出
            temp = list[i]
            # 保存下标
            index = i
            # 从后往前依次比较每个元素
            for j in range(i - 1, -1, -1):
                # 和比取出元素大的元素交换
                if list[j] > temp:
                    list[j + 1] = list[j]
                    index = j
                else:
                    break
            # 插入元素
            list[index] = temp
    return list

with open(inputfilename1,'r') as f:
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
        elif re.search('Atoms',line):
            natomstart = ncount + 2
            natomend = natomstart + natom - 1
        elif re.search('Bonds',line):
            nbondstart = ncount + 2
            nbondend = nbondstart + nbond - 1

with open(inputfilename1, 'r') as f:
    ncount = 0
    for line in f.readlines():
        ncount = ncount + 1
        if ncount >= natomstart and ncount <= natomend:
            line = line.strip().split()
            AtomID[line[0]] = [line[4], line[5], line[6]]
        if ncount >= nbondstart and ncount <= nbondend:
            line = line.strip().split()
            BondID[line[0]] = [line[2], line[3]]
print('reading file complete')

for i in range(1, int(len(AtomID)) + 1):
    temp = []
    temp.clear()
    for j in range(1, int(len(BondID)) + 1):
        for k in range(2):
            if i == int(BondID[str(j)][k]):
                temp.append(BondID[str(j)][1 - k])
    Atombond[str(i)] = temp
print('constructing dict of atombond complete')

'''
AngleIDtmp =[]
for j in range(1,len(Atombond)+1):
    G = nx.Graph(Atombond)
    p = list(nx.cycle_basis(G,str(j)))
    for i in range(len(p)):
        if len(p[i]) == 3:
            listtmp = insert_sort(p[i])
            if listtmp not in AngleIDtmp:
                AngleIDtmp.append(listtmp)
for i in range(len(AngleIDtmp)):
    AngleID[str(i+1)] = AngleIDtmp[i] 
'''

with open(inputfilename2, 'r') as f:
    ncount = 0
    for line in f.readlines():
        ncount = ncount + 1
        line = line.strip().split()
        AngleID[str(ncount)] = [line[0],line[1],line[2]]
print('constructing dict of angleid complete')

for i in range(1,len(BondID)+1):
    tmp = []
    for j in range(len(Atombond[str(BondID[str(i)][0])])):
        tmp.append(Atombond[str(BondID[str(i)][0])][j])
    for j in range(len(Atombond[str(BondID[str(i)][1])])):
        tmp.append(Atombond[str(BondID[str(i)][1])][j])
    for j in tmp:
        if j == BondID[str(i)][0] or j == BondID[str(i)][1]:
            tmp.remove(j)
    d1 = Counter(tmp)
    for key in list(d1.keys()):
        if d1[key] == 1:
            del d1[key]
            continue
    dihedrallisttmp = [BondID[str(i)][0],BondID[str(i)][1]]
    for key in list(d1.keys()):
        dihedrallisttmp.append(key)
    DihedralID[str(i)] = dihedrallisttmp
print('creating dihedralid complete')
print(DihedralID)

for key in BondID:
    p1=Point(float(AtomID[BondID[key][0]][0]),float(AtomID[BondID[key][0]][1]),float(AtomID[BondID[key][0]][2]))
    p2=Point(float(AtomID[BondID[key][1]][0]),float(AtomID[BondID[key][1]][1]),float(AtomID[BondID[key][1]][2]))
    l=Getlen(p1,p2)
    d=l.getlen()
    bondlenth[key]=[d]
print('creating bondlenth complete')

for key in AngleID:
    p1 = Point(float(AtomID[AngleID[key][0]][0]), float(AtomID[AngleID[key][0]][1]), float(AtomID[AngleID[key][0]][2]))
    p2 = Point(float(AtomID[AngleID[key][1]][0]), float(AtomID[AngleID[key][1]][1]), float(AtomID[AngleID[key][1]][2]))
    p3 = Point(float(AtomID[AngleID[key][2]][0]), float(AtomID[AngleID[key][2]][1]), float(AtomID[AngleID[key][2]][2]))
    d1 = Getlen(p1,p2)
    l1 = d1.getlen()
    d2 = Getlen(p1,p3)
    l2 = d2.getlen()
    d3 = Getlen(p2,p3)
    l3 = d3.getlen()
    a = Getarea(l1,l2,l3)
    area = a.getarea()
    anglearea[key]=[area]
print('creating anglearea complete')

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
            key + " 1001 1 " + str(AtomID[key][0]) + " " + str(AtomID[key][1]) + " " + str(AtomID[key][2]) + " 1.0 " + "\n"
        )
    f.write(
        '''
Bonds

'''
    )
    for key in BondID:
        f.write(
            key + " 1 " + str(BondID[key][0]) + " " + str(BondID[key][1]) + " " + str(bondlenth[key][0]) + "\n"
        )
    f.write(
        '''
Angles

'''
    )
    for key in AngleID:
        f.write(
            key + " 1 " + str(AngleID[key][0]) + " " + str(AngleID[key][1]) + " " + str(AngleID[key][2]) + " " + str(anglearea[key][0]) + "\n"
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

areatotal = 0
for key in anglearea:
    areatotal = areatotal + float(anglearea[key][0])
volumetotal = 0.0
gx = gy = gz = 0.0
for key in AtomID:
    gx += float(AtomID[key][0])
    gy += float(AtomID[key][1])
    gz += float(AtomID[key][2])
gx = gx / len(AtomID)
gy = gy / len(AtomID)
gz = gz / len(AtomID)
for key in AngleID:
    arr = np.array([[AtomID[AngleID[key][0]][0],AtomID[AngleID[key][1]][0],AtomID[AngleID[key][2]][0],gx],
                    [AtomID[AngleID[key][0]][1],AtomID[AngleID[key][1]][1],AtomID[AngleID[key][2]][1],gy],
                    [AtomID[AngleID[key][0]][2],AtomID[AngleID[key][1]][2],AtomID[AngleID[key][2]][2],gz],
                    [1,1,1,1]])
    arr = np.array(arr,dtype='float')
    volumetotal += 1/6*abs(np.linalg.det(arr))
print('surface area by triangle is : ',areatotal)
print('volume is : ',volumetotal)