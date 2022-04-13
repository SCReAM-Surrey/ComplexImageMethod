#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 11:42:35 2022

@author: od0014
"""

import numpy as np
import ComplexImageMethod.src.ComplexImageSimulation as sim
import ComplexImageMethod.src.Geometry as geom
from SDNPy.src.utils import visualize_room as vis
from matplotlib import cm
import matplotlib.pyplot as plt

#cuboid room dimensions
L = 4.0 
W = 5.0
H = 4.5
room = geom.Room()
room.shape = geom.Cuboid(L,W,H)


#number of wave numbers 
Nx = 50
#number of microphone positions
Ny = 50

srcPos = geom.Point(0.5,2.3,1.2)
sourceAngle = np.arctan(srcPos.y/srcPos.x)
wave_nums = np.linspace(0.01,10,Nx)

#receivers are in an arc around xz plane
r = 2;    
theta = np.linspace(0,np.pi/2,Ny)
rx = r * np.cos(theta)
ry = L/2 * np.ones([Ny,1])
rz = r * np.sin(theta)
receiverPos = list()

for i in range(Ny):
    receiverPos.append(geom.Point(rx[i],ry[i],rz[i]))

#add wall impedance
# loop over 6 walls
nWalls = room.shape.nWalls
# constant impedance over all frequencies


for k in range(nWalls):
    room.wallImpedance.append([0.95+1j*0.5 for i in range(Nx)])

#####################################
# draw the setup
# plot to visualize room - optional
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
vis.plot_room(ax,
              [room.shape.x / 2, room.shape.y / 2, room.shape.z / 2],
              room.shape.x, room.shape.y, room.shape.z)
vis.plot_point(ax, srcPos, 'source')
for k in range(Ny):
    vis.plot_point(ax, receiverPos[k], 'mic')
ax.view_init(45,110)
plt.savefig('../figures/test1_setup.png')
plt.show()


######################################
#run the simulation
order = 3
csim = sim.ComplexImageSimulation(room, srcPos, receiverPos, wave_nums, order)
ref_pressure, total_pressure = csim.run()


######################################
# surface plot
[K, Theta] = np.meshgrid(wave_nums,theta)
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
surf = ax.plot_surface(K, Theta/np.pi, 20*np.log10(np.abs(total_pressure)), cmap=cm.jet, linewidth=0.1)
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('Wave number')
ax.set_ylabel('Angle of receiver from origin')
ax.set_zticks([])
ax.view_init(270,-90)
plt.savefig('../figures/test1_order=' + str(order) + '_surf.png')
plt.show()


#polar plot
idx = [1,5,10]
fig = plt.figure()
for j in range(len(idx)):
    polar_p = total_pressure[:,idx[j]];
    ax = fig.add_subplot(len(idx),1,j+1, polar=True)
    ax.plot(theta, np.abs(polar_p));
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    ax.set_title('k = ' + str(round(wave_nums[idx[j]],3)), fontdict={'fontsize': 8, 'fontweight': 'light'})
    
plt.savefig('../figures/test1_order=' + str(order) + '_polar.png')
plt.show()


    


