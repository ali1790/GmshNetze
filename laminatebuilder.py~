#!/usr/bin/python3.8
import glob
import math
import matplotlib.pyplot as plt
from noise import pnoise1, pnoise3
import numpy as np
import os
import random
import sys
#Perlin noise is a random signal that is correlated over short time scales
# The number of octaves determines the time/length scales at which features in the signal are present
# The persistence determines the amplitude of the features
class Geowriter():
    def __init__(self):
        self.pointid = 0
        self.lineid = 0
        self.loopid = 0
        self.surfaceid = 0
        self.file = ''
        self.initgmsh()

    def initgmsh(self):
        self.gmsh_link = os.path.islink('gmsh')

        if not self.gmsh_link:
            print ('gmsh ist nicht verlinkt')
            from shutil import which
            if which('gmsh') is not None:
                print ('gmsh ist installiert')
                self.gmsh_inst = True
            else:
                sys.exit('gmsh ist weder installiert noch verlinkt')

    def mesh(self):
        if self.gmsh_link:
            os.system('./gmsh -2 %s' % (self.file))
        elif self.gmsh_inst:
            os.system('gmsh -2 %s' % (self.file))

    def convert(self):
        os.system('dolfin-convert %s %s' % (self.file.replace('.geo', '.msh'), self.file.replace('.geo', '.xml')))
    def reverselines(self, linelist):
        return [-l for l in linelist[::-1]]

    def writecircle(self, p1, pc, p2):
        with open(self.file, 'a') as o:
            self.lineid +=1
            o.write('Circle(%i) = {%i, %i, %i};\n' % (self.lineid, p1, pc, p2))
        return self.lineid

    def writeline(self, p1, p2):
        with open(self.file, 'a') as o:
            self.lineid +=1
            o.write('Line(%i) = {%i, %i};\n' % (self.lineid, p1, p2) )
        return self.lineid

    def writelineloop(self, linelist, s = 'u'):
        with open(self.file, 'a') as o:
            self.loopid += 1
            o.write('Line Loop(%i) = {%s};\n' % (self.loopid, ','.join([str(s) for s in linelist])))
        self.writeplanesurface(self.loopid)
        if s == 's':
            self.writetransfinitesurface(self.loopid)
        return self.loopid

    def writepoint(self, v, mesh):
        with open(self.file, 'a') as o:
            self.pointid +=1
            o.write('Point(%i) = {%f, %f, 0.0, %f};\n' % (self.pointid, v[0], v[1], mesh))
        return self.pointid

    def writeplanesurface(self, loopid):
        with open(self.file, 'a') as o:
            self.surfaceid += 1
            o.write('Plane Surface(%i) = {%i};\n' % (self.surfaceid, loopid))
        return self.surfaceid

    def writetransfinitesurface(self, loopid):
        with open(self.file, 'a') as o:
            o.write('Transfinite Surface {%i};\n' % (loopid))

    @staticmethod
    def createperlin(npoints, l, noct, pers, hmax):
        ran = int(random.random() * l)
        x = np.linspace(ran, ran + l, num = npoints)
        signal = np.zeros(npoints)
        for ctr, t in enumerate(x):
            signal[ctr] = pnoise1(t, octaves = noct, persistence = pers)
        scalefactor = hmax / max([max(signal), abs(min(signal))])
        datastr = 'npoints = %i, l = %f, noct = %i, pers = %f, hmax = %f' % (npoints, l, noct, pers, hmax)
        return [[[t - min(x), s * scalefactor] for t, s in zip(x, signal)], datastr]

class ClassicLaminate(Geowriter):
    def __init__(self, l, h1, h2, dh, cm, fm, filename = 'perlinlaminate.geo'):
        super().__init__()
        self.file = filename
        print (self.file)
        
        with open(filename, 'w') as o:
            o.write('//foo\n')


        p1 = self.writepoint([0.0, 0.0], cm)
        p2 = self.writepoint([l, 0.0], cm)
        p3 = self.writepoint([0.0, h1], fm)
        p4 = self.writepoint([l, h1], fm)

        l1 = self.writeline(p1, p2)
        l2 = self.writeline(p1, p3)
        l3 = self.writeline(p2, p4)
        l4 = self.writeline(p3, p4)

        ll1 = self.writelineloop([l1, l3, -l4, -l2])

        p5 = self.writepoint([0.0, h1 +dh], fm)
        p6 = self.writepoint([l, h1 + dh], fm)

        l5 = self.writeline(p3, p5)
        l6 = self.writeline(p4, p6)
        l7 = self.writeline(p5, p6)

        ll2 = self.writelineloop([l4, l6, -l7, -l5])

        p7 = self.writepoint([0.0, h1 + h2 -dh], fm)
        p8 = self.writepoint([l, h1 + h2 -dh], fm)
        
        l8 = self.writeline(p5, p7)
        l9 = self.writeline(p6, p8)
        l10 = self.writeline(p7, p8)
        
        ll3 = self.writelineloop([l7, l9, -l10, -l8])
    
        p9  = self.writepoint([0.0, h1 + h2], fm)
        p10 = self.writepoint([l, h1 + h2], fm)

        l11 = self.writeline(p7, p9)
        l12 = self.writeline(p8, p10)
        l13 = self.writeline(p9, p10)

        ll4 = self.writelineloop([l10, l12, -l13, -l11])

        p11 = self.writepoint([0.0, 2.0 * h1 + h2], cm)
        p12 = self.writepoint([l, 2.0 * h1 + h2], cm)

        l14 = self.writeline(p9, p11)
        l15 = self.writeline(p10, p12)
        l16 = self.writeline(p11, p12)

        ll5 = self.writelineloop([l13, l15, -l16, -l14])

        with open(filename, 'a') as o:
            o.write('Physical Surface(1) = {1};')
            o.write('Physical Surface(2) = {2};')
            o.write('Physical Surface(3) = {3};')
            o.write('Physical Surface(4) = {4};')
            o.write('Physical Surface(5) = {5};')

class AngledLaminate(Geowriter):
    def __init__(self, l, h1, h2, a1, a2, dh, cm, fm, filename = 'perlinlaminate.geo'):
        super().__init__()
        self.file = filename
        
        alpha1 = a1 * (np.pi / 180.0)
        alpha2 = a2 * (np.pi / 180.0)

        dh1 = np.tan(alpha1) * l
        dh2 = np.tan(alpha2) * l

        with open(filename, 'w') as o:
            o.write('//foo\n')


        p1 = self.writepoint([0.0, 0.0], cm)
        p2 = self.writepoint([l, 0.0], cm)
        p3 = self.writepoint([0.0, h1], fm)
        p4 = self.writepoint([l, h1 + dh1], fm)

        l1 = self.writeline(p1, p2)
        l2 = self.writeline(p1, p3)
        l3 = self.writeline(p2, p4)
        l4 = self.writeline(p3, p4)

        ll1 = self.writelineloop([l1, l3, -l4, -l2])

        p5 = self.writepoint([0.0, h1 +dh], fm)
        p6 = self.writepoint([l, h1 + dh + dh1], fm)

        l5 = self.writeline(p3, p5)
        l6 = self.writeline(p4, p6)
        l7 = self.writeline(p5, p6)

        ll2 = self.writelineloop([l4, l6, -l7, -l5])

        p7 = self.writepoint([0.0, h1 + h2 -dh], fm)
        p8 = self.writepoint([l, h1 + h2 -dh + dh2], fm)
        
        l8 = self.writeline(p5, p7)
        l9 = self.writeline(p6, p8)
        l10 = self.writeline(p7, p8)
        
        ll3 = self.writelineloop([l7, l9, -l10, -l8])
    
        p9  = self.writepoint([0.0, h1 + h2], fm)
        p10 = self.writepoint([l, h1 + h2 + dh2], fm)

        l11 = self.writeline(p7, p9)
        l12 = self.writeline(p8, p10)
        l13 = self.writeline(p9, p10)

        ll4 = self.writelineloop([l10, l12, -l13, -l11])

        p11 = self.writepoint([0.0, 2.0 * h1 + h2], cm)
        p12 = self.writepoint([l, 2.0 * h1 + h2], cm)

        l14 = self.writeline(p9, p11)
        l15 = self.writeline(p10, p12)
        l16 = self.writeline(p11, p12)

        ll5 = self.writelineloop([l13, l15, -l16, -l14])

        with open(filename, 'a') as o:
            o.write('Physical Surface(1) = {1};')
            o.write('Physical Surface(2) = {2};')
            o.write('Physical Surface(3) = {3};')
            o.write('Physical Surface(4) = {4};')
            o.write('Physical Surface(5) = {5};')


class PerlinLaminate(Geowriter):
    def __init__(self, S1 ,S2 ,h1 ,h2 ,cm ,fm ,dh ,filename = 'perlinlaminate.geo'):
        super().__init__()
        self.s1 = [[v[0], v[1] + h1] for v in S1[0]] # Signal der unteren Grenzflaeche
        self.s1a = [[v[0], v[1] + h1 + dh] for v in S1[0]] # Signal der unteren Grenzflaeche
        self.s2 = [[v[0], v[1] + h2 + h1 - dh] for v in S2[0]] # Signal der oberen Grenzflaeche
        self.s2a = [[v[0], v[1] + h2 + h1] for v in S2[0]] # Signal der oberen Grenzflaeche
        self.h1 = h1 # Hoehe der unteren und oberen Schicht
        self.h2 = h2 # Hoehe der mittleren Schicht
        self.cmesh = cm # Grobes Netz
        self.fmesh = fm # Feines Netz
        self.file = filename
        with open(self.file, 'w') as o: 
            o.write('//Laminat mit rauen Grenzflaechen\n')
            o.write('//L = %f\n' % (self.s1[0][-1] - self.s1[0][0]))
            o.write('//h1 = %f, h2 = %f \n' % (self.h1, self.h2))
            o.write('//Signalparameter unten:\n')
            o.write('//%s\n' % (S1[1]))
            o.write('//Signalparameter oben:\n')
            o.write('//%s\n' % (S2[1]))

        p1 = self.writepoint([0.0, 0.0], self.cmesh)
        p2 = self.writepoint([self.s1[-1][0], 0.0], self.cmesh)
        interfacepoints1 = []
        interfacelines1 = []
        for i, v in enumerate(self.s1):
            interfacepoints1.append(self.writepoint(v, self.fmesh))
            if i>0:
               interfacelines1.append(self.writeline(interfacepoints1[-2], interfacepoints1[-1]))

        l1 = self.writeline(p1, interfacepoints1[0])
        l2 = self.writeline(p1, p2)
        l3 = self.writeline(p2, interfacepoints1[-1])
        ll1 = self.writelineloop([-l1, l2, l3] + self.reverselines(interfacelines1))

        interfacepoints1a = []
        interfacelines1a = []

        for i, v in enumerate(self.s1a):
            interfacepoints1a.append(self.writepoint(v, self.fmesh))
            if i>0:
               interfacelines1a.append(self.writeline(interfacepoints1a[-2], interfacepoints1a[-1]))
        l1a = self.writeline(interfacepoints1[0], interfacepoints1a[0])
        l2a = self.writeline(interfacepoints1[-1], interfacepoints1a[-1])
        ll1a = self.writelineloop(interfacelines1 + [l2a] + self.reverselines(interfacelines1a) + [-l1a])

        interfacepoints2 = []
        interfacelines2 = []
        for i, v in enumerate(self.s2):
            interfacepoints2.append(self.writepoint(v, self.fmesh))
            if i>0:
               interfacelines2.append(self.writeline(interfacepoints2[-2], interfacepoints2[-1]))

        l1 = self.writeline(interfacepoints1a[0], interfacepoints2[0])
        l2 = self.writeline(interfacepoints1a[-1], interfacepoints2[-1])
        ll2 = self.writelineloop(interfacelines1a + [l2] + self.reverselines(interfacelines2) + [-l1])

        interfacepoints2a = []
        interfacelines2a = []
        for i, v in enumerate(self.s2a):
            interfacepoints2a.append(self.writepoint(v, self.fmesh))
            if i>0:
               interfacelines2a.append(self.writeline(interfacepoints2a[-2], interfacepoints2a[-1]))

        l1 = self.writeline(interfacepoints2[0], interfacepoints2a[0])
        l2 = self.writeline(interfacepoints2[-1], interfacepoints2a[-1])
        ll2a = self.writelineloop(interfacelines2 + [l2] + self.reverselines(interfacelines2a) + [-l1])

        p1 = self.writepoint([0.0, 2.0 * h1 + h2], self.cmesh)
        p2 = self.writepoint([self.s1[-1][0], 2.0 * h1 + h2], self.cmesh)
        l1 = self.writeline(interfacepoints2a[0], p1)
        l2 = self.writeline(p1, p2)
        l3 = self.writeline(interfacepoints2a[-1], p2)
        ll3 = self.writelineloop(interfacelines2a + [l3, -l2, -l1])

class CurvedLaminate(Geowriter):
    def __init__(self ,l , h1 ,h2 ,dh1 ,dh2 ,cm ,fm ,dh ,filename = 'curvedlaminate.geo'):
        super().__init__()

        self.file = filename

        with open(filename, 'w') as o:
            o.write('//Laenge = %f\n' % (l))
            o.write('//h1 = %f\n' % (h1))
            o.write('//h2 = %f\n' % (h2))
            o.write('//dh1 = %f\n' % (dh1))
            o.write('//dh2 = %f\n' % (dh2))
            o.write('//r1 = %f\n' % ((4.0 * (dh1 ** 2.0) + l ** 2.0) / (8.0 * dh1)))
            o.write('//r2 = %f\n' % ((4.0 * (dh2 ** 2.0) + l ** 2.0) / (8.0 * dh2)))
            o.write('//foo\n')
        # Daten des Ersten Kreises
        if dh1 > 0.0:
            print ('Test')
            v_center1 = [0.5 * l, h1 - (((l ** 2.0) - 4.0 * (dh1 ** 2.0)) / (8.0 * abs(dh1)))]
            v_center1A = [0.5 * l, h1 +dh - (((l ** 2.0) - 4.0 * (dh1 ** 2.0)) / (8.0 * abs(dh1)))]
            alpha1 = 4.0 * np.arctan2(2.0 * dh1 , l)
        else:
            print ('Test2')
            v_center1 = [0.5 * l, h1 + (((l ** 2.0) - 4.0 * (dh1 ** 2.0)) / (8.0 * abs(dh1)))]
            v_center1A = [0.5 * l, h1 +dh + (((l ** 2.0) - 4.0 * (dh1 ** 2.0)) / (8.0 * abs(dh1)))]
            alpha1 = 4.0 * np.arctan2(2.0 * dh1 , l)

        if dh2 > 0:
            v_center2 = [0.5 * l, h1 + h2 - (((l ** 2.0) - 4.0 * (dh2 ** 2.0)) / (8.0 * abs(dh2)))]
            v_center2A = [0.5 * l, h1 + h2 + dh - (((l ** 2.0) - 4.0 * (dh2 ** 2.0)) / (8.0 * abs(dh2)))]
            alpha2 = 4.0 * np.arctan2(2.0 * dh2 , l)
        else:
            v_center2 = [0.5 * l, h1 + h2 + (((l ** 2.0) - 4.0 * (dh2 ** 2.0)) / (8.0 * abs(dh2)))]
            v_center2A = [0.5 * l, h1 + h2 + dh + (((l ** 2.0) - 4.0 * (dh2 ** 2.0)) / (8.0 * abs(dh2)))]
            alpha2 = 4.0 * np.arctan2(2.0 * dh2 , l)

        p0 = self.writepoint(v_center1, fm)
        p01 = self.writepoint(v_center1A, fm)
        p1 = self.writepoint([0.0, 0.0], cm)
        p2 = self.writepoint([l, 0.0], cm)
        p3 = self.writepoint([0.0, h1], fm)
        p4 = self.writepoint([l, h1], fm)
        p5 = self.writepoint([0.5 * l, h1 + dh1], fm)

        l1 = self.writeline(p1, p2)
        l2 = self.writeline(p1, p3)
        l3 = self.writeline(p2, p4)
        l4 = self.writecircle(p3, p0, p5)
        l5 = self.writecircle(p5, p0, p4)

        ll1 = self.writelineloop([l1, l3, -l5, -l4, -l2])
        
        p6 = self.writepoint([0.0, h1 + dh], fm)
        p7 = self.writepoint([l, h1 + dh], fm)
        p8 = self.writepoint([0.5 * l, h1 + dh1 + dh], fm)
        
        l6 = self.writeline(p3, p6)
        l7 = self.writeline(p4, p7)
        l8 = self.writecircle(p6, p01, p8)
        l9 = self.writecircle(p8, p01, p7)

        ll2 = self.writelineloop([l4, l5, l7, -l9, -l8, -l6])

        p0 = self.writepoint(v_center2A, fm)
        p9 = self.writepoint([0.0, h1 + h2 - dh], fm)
        p10 = self.writepoint([l, h1 + h2 - dh], fm)
        p11 = self.writepoint([0.5 * l, h1 + h2 - dh + dh2], fm)

        l10 = self.writeline(p6, p9)
        l11 = self.writeline(p7, p10)
        l12 = self.writecircle(p9, p0, p11)
        l13 = self.writecircle(p11, p0, p10)
        
        ll3 = self.writelineloop([l8, l9, l11, -l13, -l12, -l10])

        p12 = self.writepoint([0.0, h1 + h2], fm)
        p13 = self.writepoint([l , h1 + h2], fm)
        p14 = self.writepoint([0.5 * l, h1 + h2 + dh2], fm)

        l14 = self.writeline(p9, p12)
        l15 = self.writeline(p10, p13)
        l16 = self.writecircle(p12, p0, p14)
        l17 = self.writecircle(p14, p0, p13)

        ll4 = self.writelineloop([l12, l13, l15, -l17, -l16, -l14])

        p15 = self.writepoint([0.0, 2.0 * h1 + h2], cm)
        p16 = self.writepoint([l, 2.0 * h1 + h2], cm)

        l18 = self.writeline(p12, p15)
        l19 = self.writeline(p13, p16)
        l20 = self.writeline(p15, p16)

        ll4 = self.writelineloop([l16, l17, l19, -l20, -l18])

def isfloat(val):
    try:
        float(val)
        return True
    except ValueError:
        return False

def readoptions():
    tmp = {}
    with open('laminate.dat', 'r') as i: lines = i.readlines()
    for l in lines:
        if '#' not in l and len(l) > 1 :
            key = l.rstrip().replace(' ','').split('=')[0]
            val = l.rstrip().replace(' ','').split('=')[1]
            if isfloat(val):
                tmp[key] = float(val)
            else:
                tmp[key] = str(val)

    return tmp

opt = readoptions()
if opt['Typ'] == 'Rect':
    print ('Typ: Rect')
    name = opt['name']
    if '.geo' not in name:
        name = name +'.geo'
    print ('Dateiname '+ name + '.geo')
    l = opt['Laenge']
    print ('Laenge: %f' % (l))
    h1 = opt['Hoehe1']
    print ('Hoehe der aeusseren Schichten: %f' % (h1))
    h2 = opt['Hoehe2']
    print ('Hoehe der inneren Schicht: %f' % (h1))
    fm = opt['finemesh']
    print ('Feine Vernetzung: %f' % (fm))
    cm = opt['coarsemesh']
    print ('Grobe Vernetzung: %f' % (cm))
    dh = opt['h_interface']
    print ('Hoehe der Grenzfläche: %f' % (dh))
    
    test = ClassicLaminate(l, h1, h2, dh, cm, fm, filename = name)

elif opt['Typ'] == 'Angled':
    print ('Typ: Angled')
    name = opt['name']
    if '.geo' not in name:
        name = name +'.geo'
    print ('Dateiname '+ name + '.geo')
    l = opt['Laenge']
    print ('Laenge: %f' % (l))
    h1 = opt['Hoehe1']
    print ('Hoehe der aeusseren Schichten: %f' % (h1))
    h2 = opt['Hoehe2']
    print ('Hoehe der inneren Schicht: %f' % (h1))
    alpha1 = opt['alpha1']
    print ('Winkel der unteren Schicht')
    alpha2 = opt['alpha2']
    print ('Winkel der oberen Schicht')
    fm = opt['finemesh']
    print ('Feine Vernetzung: %f' % (fm))
    cm = opt['coarsemesh']
    print ('Grobe Vernetzung: %f' % (cm))
    dh = opt['h_interface']
    print ('Hoehe der Grenzfläche: %f' % (dh))
    
    test = AngledLaminate(l, h1, h2, alpha1, alpha2, dh, cm, fm, filename = name)
elif opt['Typ'] == 'Curved':
    name = opt['name']
    print ('Dateiname '+ name + '.geo')
    print ('Typ: Curved')
    l = opt['Laenge']
    print ('Laenge: %f' % (l))
    h1 = opt['Hoehe1']
    print ('Hoehe der aeusseren Schichten: %f' % (h1))
    h2 = opt['Hoehe2']
    print ('Hoehe der inneren Schicht: %f' % (h1))
    dh1 = opt['dh1']
    print ('Maximale Abweichung der unteren Schicht : %f' % (dh1))
    dh2 = opt['dh2']
    print ('Maximale Abweichung der oberen Schicht : %f' % (dh2))
    fm = opt['finemesh']
    print ('Feine Vernetzung: %f' % (fm))
    cm = opt['coarsemesh']
    print ('Grobe Vernetzung: %f' % (cm))
    dh = opt['h_interface']
    print ('Hoehe der Grenzfläche: %f' % (dh))
    test = CurvedLaminate(l, h1, h2, dh1, dh2, cm, fm, dh, filename = name)
