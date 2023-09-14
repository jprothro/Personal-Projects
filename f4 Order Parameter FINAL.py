#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 09:47:33 2023

@author: jadenprothro

This program is able to calculate the average f4 order parameter of a given trajectory for each frame and will also
present a graph showing the general trend of the f4 order parameter to see if the group of atoms is forming a
clathrate hydrate structure. You will need a .gro file or trajectory to run this program.
"""

import MDAnalysis
import numpy as np
import math
import matplotlib.pyplot as plt

u = MDAnalysis.Universe("in.gro")
b = MDAnalysis.Universe("in.xtc")
vdwradii = {'H':2.1, 'DUMMY':1.1, 'O':1.5}
tsangles = []
f4averages = []
frames = []

for ts in b.trajectory[:2000:500]:
    angles = []
    f4List = []
    sol = u.select_atoms("resname SOL")
    u.atoms.positions = ts.positions
    residueNum = sol.residues
    
    #goes through every single residue and calculates the angle between the planes of any bonded molecules of each residue
    for k in range(len(residueNum)): # goes through every residue in the universe
        name1 = "byres (around 3 (resnum " + str(k+1) + " and name OW)) and name OW"
        name2 = "resnum " + str(k+1)
        nearbyOxygen = sol.select_atoms(name1) #Finds all of the residues of Oxygen within 2.8 angstroms of residue 1
        mainRes = sol.select_atoms(name2)
        mainAtom = mainRes.select_atoms("name OW")
        positions1 = [ts.positions[mainAtom.ids[0]-1], ts.positions[mainAtom.ids[0]], ts.positions[mainAtom.ids[0]+1] ] # getting the x,y,z positions of the oxygen, hydrogen1, and hydrogen2 respectively in the main residue from the trajectory file
        i = 0
        for residues in nearbyOxygen.residues:

            angles2 = []
            secondMainAtom = residues.atoms.select_atoms("name OW")
            positions2 = [ts.positions[secondMainAtom.ids[0]-1], ts.positions[secondMainAtom.ids[0]], ts.positions[secondMainAtom.ids[0]+1] ]
            
            # finding the outermost hydrogens
            # distance equation: distance = math.sqrt(((positions[0][0]-x))**2 + ((positions[0][1]-y))**2 + ((positions[0][2]-z))**2)
            distance1H1 = math.sqrt(((positions1[1][0]-positions2[0][0]))**2 + ((positions1[1][1]-positions2[0][1]))**2 + ((positions1[1][2]-positions2[0][2]))**2)
            distance1H2 = math.sqrt(((positions1[2][0]-positions2[0][0]))**2 + ((positions1[2][1]-positions2[0][1]))**2 + ((positions1[2][2]-positions2[0][2]))**2)
            distance2H1 = math.sqrt(((positions2[1][0]-positions1[0][0]))**2 + ((positions2[1][1]-positions1[0][1]))**2 + ((positions2[1][2]-positions1[0][2]))**2)
            distance2H2 = math.sqrt(((positions2[2][0]-positions1[0][0]))**2 + ((positions2[2][1]-positions1[0][1]))**2 + ((positions2[2][2]-positions1[0][2]))**2)
            
            if distance1H1 > distance1H2: # picking the hydrogen furthest outside
                vector1 = [(positions1[0][0]-positions1[1][0]), (positions1[0][1]-positions1[1][1]), (positions1[0][2]-positions1[1][2])] #Main oxygen to hydrogen 1
            else:
                vector1 = [(positions1[0][0]-positions1[2][0]), (positions1[0][1]-positions1[2][1]), (positions1[0][2]-positions1[2][2])] #Main oxygen to hydrogen 1
            
            vector2 = [(positions1[0][0]-positions2[0][0]), (positions1[0][1]-positions2[0][1]), (positions1[0][2]-positions2[0][2])] #Main oxygen to second main oxygen
            
            if distance2H1 > distance2H2: # picking the hydrogen furthest outside
                vector3 = [(positions2[0][0]-positions2[1][0]), (positions2[0][1]-positions2[1][1]), (positions2[0][2]-positions2[1][2])] #second main oxygen to hydrogen 1
            else:
                vector3 = [(positions2[0][0]-positions2[2][0]), (positions2[0][1]-positions2[2][1]), (positions2[0][2]-positions2[2][2])] #second main oxygen to hydrogen 1
                
            vector4 = [(positions2[0][0]-positions1[0][0]), (positions2[0][1]-positions1[0][1]), (positions2[0][2]-positions1[0][2])] #second main oxygen to main oxygen
            
            #crossing the vectors of the respective molecule to find the normal vector of each molecule's plane
            n1 = np.cross(vector1, vector2)
            n2 = np.cross(vector3, vector4)
            
            #finding the magnitude of each molecules normal vector
            magnitude1 = math.sqrt(n1[0]**2 + n1[1]**2 + n1[2]**2)
            magnitude2 = math.sqrt(n2[0]**2 + n2[1]**2 + n2[2]**2)
            
            #calculating the angle between the two planes in radians
            planeAngle = np.arccos(np.dot(n2, n1)/(magnitude1*magnitude2))
            
            f4 = math.cos(3*(planeAngle-3.14159263))
            angles2.append(k)
            angles2.append(residues.resnum)
            angles2.append(i)
            angles2.append(planeAngle)
            angles2.append(f4)
            f4List.append(f4)
            i += 1
            angles.append(angles2)
    # weird system of adding to avoid null values in the list
    f4Sum = 0;
    for f4 in f4List:
        if math.isnan(f4) == False:
            f4Sum += f4
        else:
            f4Sum +=0
    f4averages.append((f4Sum)/len(f4List))
    frames.append(ts.frame)
    lenf4 = len(f4List)
    tsangles.append((f4Sum)/len(f4List))
    tsangles.append(angles)       
    print(f"avg f4 for frame {ts.frame} {(f4Sum)/len(f4List):.8f}")

plt.plot(frames,f4averages)
plt.xlabel("Frame")
plt.ylabel("f4 Order Parameter Average")
plt.show()
