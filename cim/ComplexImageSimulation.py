#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 11:36:46 2022

@author: od0014
"""

import numpy as np
import cim.SimpleImageSource as simp_img


class ComplexImageSimulation:
    def __init__(self, room, source, mic_array, wave_num, order, *args):
        self.room = room
        self.source = source
        self.mic_array = mic_array
        self.wave_num = wave_num
        self.maxOrder = order  # reflection order
        if len(args) == 0:
            self.finiteWall = False
        else:
            self.finiteWall = args[0]

    def run(self):

        walls = self.room.shape.setWallPosition(self.room.wallImpedance)
        # total_IM = 0

        # I am trying to build a dictionary with tuple as keys
        wallNames = [*walls]
        IDs = list(tuple())
        for j in range(self.room.shape.nWalls):
            for i in range(1, self.maxOrder + 2):
                IDs.append((wallNames[j], i))

        self.imageSrc = {(wallID, order): [] for (wallID, order) in IDs}

        reflectedPressure = np.zeros(
            (len(self.mic_array), len(self.wave_num)), dtype=complex
        )

        for i in range(1, self.maxOrder + 1):

            imagesByOrder = []

            for wallID in wallNames:

                wall = walls[wallID]

                if i == 1:

                    S0 = simp_img.ImageSource(
                        self.source,
                        1.0,
                        wall.plane.getPointReflection(self.source),
                        wall,
                        i,
                    )
                    # find receiver distance and angle from image source
                    S0.getRelativeReceiverParameters(self.mic_array)
                    # calculate image strength
                    S0.calculateImageStrength(wall, self.wave_num)

                    self.imageSrc[wallID, i].append(S0)

                for S in self.imageSrc[wallID, i]:
                    # calculate reflected pressure due to this image source

                    # add reflcted pressure only if that image source does not exist
                    existFlag = [S.pos.equals(im.pos) for im in imagesByOrder]
                    if True in existFlag:
                        continue
                    else:
                        [k, r] = np.meshgrid(self.wave_num, S.r)
                        kr = np.multiply(k, r)
                        reflectedPressure += np.divide(np.exp(-1j * kr), r) * S.strength

                        imagesByOrder.append(S)
                        # total_IM += 1

                        if i < self.maxOrder:
                            for otherWallID in wallNames:

                                if otherWallID == wallID:
                                    continue
                                else:
                                    otherWall = walls[otherWallID]

                                    # check if current image source is behind wall or not
                                    if not otherWall.isPointBehindWall(S.pos):

                                        Sm = simp_img.ImageSource(
                                            S.pos,
                                            S.strength,
                                            otherWall.plane.getPointReflection(S.pos),
                                            otherWall,
                                            i + 1,
                                        )

                                        Sm.getRelativeReceiverParameters(self.mic_array)

                                        Sm.calculateImageStrength(
                                            otherWall, self.wave_num
                                        )

                                        self.imageSrc[otherWallID, i + 1].append(Sm)

            print("Reflected pressure calculated for reflection order ", str(i))

        # calculate direct pressure
        dis = np.zeros(len(self.mic_array))
        for i in range(len(self.mic_array)):
            dis[i] = self.source.getDistance(self.mic_array[i])

        [k, r0] = np.meshgrid(self.wave_num, dis)
        kr0 = np.multiply(k, r0)
        directPressure = np.divide(np.exp(-1j * kr0), r0)

        # total pressure is the sum of direct and reflceted pressures
        totalPressure = directPressure + reflectedPressure

        # print('Total number of image sources is ' , total_IM)
        return reflectedPressure, totalPressure
