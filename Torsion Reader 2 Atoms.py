#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 20:16:20 2023

@author: jadenprothro

This program can calculate the bond angle torsion between 2 water molecules given the x,y,z coordinates of the atoms
in a space or the positions relative to each other.
"""

import MDAnalysis
import numpy as np
import math

u = MDAnalysis.Universe("in.gro")
vdwradii = {'H':2.1, 'DUMMY':1.1, 'O':1.5}

sol = u.select_atoms("resname SOL")
# original iteration for two singular oxygen atoms
#-132.383 degrees should be the answer
#two water molecules that are bonded
OW1 = sol.select_atoms("resnum 1") 
OW7449 = sol.select_atoms("resnum 1863")

# getting the x,y,z positions of all the atoms in each molecule
positions1 = OW1.positions
positions2 = OW7449.positions

# finding the outermost hydrogens
# distance equation: distance = math.sqrt(((positions[0][0]-x))**2 + ((positions[0][1]-y))**2 + ((positions[0][2]-z))**2)
distance1H1 = math.sqrt(((positions1[1][0]-positions2[0][0]))**2 + ((positions1[1][1]-positions2[0][1]))**2 + ((positions1[1][2]-positions2[0][2]))**2)
distance1H2 = math.sqrt(((positions1[2][0]-positions2[0][0]))**2 + ((positions1[2][1]-positions2[0][1]))**2 + ((positions1[2][2]-positions2[0][2]))**2)
distance2H1 = math.sqrt(((positions2[1][0]-positions1[0][0]))**2 + ((positions2[1][1]-positions1[0][1]))**2 + ((positions2[1][2]-positions1[0][2]))**2)
distance2H2 = math.sqrt(((positions2[2][0]-positions1[0][0]))**2 + ((positions2[2][1]-positions1[0][1]))**2 + ((positions2[2][2]-positions1[0][2]))**2)


#finding the vectors from the Oxygen molecule to each hydrogen molecule for both molecules
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

#O1H1 and O2H1: 2.699
#O1H2 and O2H1: .7899
#O1H2 and O2H2: 2.7399
#*O1H1 and O2H2: .8311* WINNER



