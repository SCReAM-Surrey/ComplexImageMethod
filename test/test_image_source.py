#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 14:40:26 2022

@author: od0014
"""

import numpy as np
import ComplexImageMethod.src.ImageSource as img
import ComplexImageMethod.src.Geometry as geom

#cuboid room dimensions
L = 4.0 
W = 5.0
H = 4.5
room = geom.Room()
room.shape = geom.Cuboid(L,W,H)


#number of wave numbers 
Nx = 50
#number of microphone positions
Ny = 100

srcPos = geom.Point(0.5,2.3,1.2)
wave_nums = np.linspace(0.01,10,Nx)


h = 0.5 
receiverPos = list()
receiverPos.append(geom.Point(L/2, W/2, h))
receiverPos.append(geom.Point(4, 3, H/2))


#add wall impedance
# loop over 6 walls
nWalls = room.shape.nWalls
# constant impedance over all frequencies
for k in range(nWalls):
    room.wallImpedance.append([0.95+1j*0.5 for i in range(Nx)])


#add walls to the room
order = 1
finiteWall = False
wallList = room.shape.setWallPosition(order, room.wallImpedance)
srcStrength = {
                "floor": 1.0,
                "ceiling": 1.0,
                "left": 1.0,
                "right": 1.0,
                "front": 1.0,
                "back": 1.0}

#the following is supposed to work for first-order reflections only
for wallName in wallList:
    
    print(wallName)
    curWall = wallList[wallName][0] 

    IS = img.ImageSource(srcPos, srcStrength, curWall, order, finiteWall)
    angleLims = IS.calculateAngleLimitsWithOtherWalls(curWall)
    print(angleLims)
    
    IS.getRelativeReceiverParameters(receiverPos)
    IS.calculateImageStrength(curWall, wallName, wave_nums)
    print(IS.strength[wallName])