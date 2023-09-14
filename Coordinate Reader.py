#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 09:52:15 2022

@author: jadenprothro

This program can take a still .gro file and list out the coordinates of all the atoms or molecules of a certain type 
in an x,y,z coordinate format.
"""


import MDAnalysis
import math

u = MDAnalysis.Universe("init-zero.gro")
atomsO = u.select_atoms('name OW')
positions =  atomsO.positions
print(positions[0][0])
distances = []
i = 0

while i < len(positions):
    distances.append([])
    x = positions[i][0]
    y = positions[i][1]
    z = positions[i][2]
    for molecule in range(len(positions)):
        distance = round(math.sqrt(((positions[molecule][0]-x))**2 + ((positions[molecule][1]-y))**2 + ((positions[molecule][2]-z))**2), 3)
        distances[i].append(distance)
    i += 1 

print(distances[0][1])


'''
with open("distancesO.csv", "w") as file:
    fileC = csv.writer(file)
    for molecule in range(len(positions)):
        fileC.write("OW")
        fileC.write(f"{molecule}")
        x = positions[molecule,0]
        y = positions[molecule,1]
        z = positions[molecule,2]
        
        for line in positions:
            distance = math.sqrt((line[0]-x)**2 + (line[1]-y)**2 + (line[2]-z)**2)
            file.write(f"{distance:.0f}")
        
    fileC.writerow([])
    '''

#names_positions = [[at] + list(pos) for at, pos in zip(atomsO.names, atomsO.positions)]
#names_positions =  atomsO.positions

# show the first three entries
#print(names_positions[:3])
    
#atomsO.write('posO.gro')

'''
with MDAnalysis.Writer("protein.xtc", protein.n_atoms) as W:
    for ts in u.trajectory:
        W.write(protein)
'''