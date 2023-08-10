#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 19:20:15 2022

@author: od0014
"""
import numpy as np
import mpmath as mp
from scipy.special import erfc


def myGammainc(c, B, t):
    ypts = []
    B_list = B.tolist()
    c_list = c.to_list()
    t_list = t.to_list()
    for i in range(len(c_list)):
        for k in range(len(c_list[i])):
            ypts.append(
                mp.gammainc(0.5, c_list[i][k] * (t_list[i][k] - 1j * B_list[i][k]))
            )

    return np.reshape(np.array(ypts, dtype=complex), np.shape(c))


class ImageSource:

    """
    Class that keeps track of the image sources and their orders
    """

    def __init__(self, sourcePos, sourceStrength, pos, wall, order, img_type : str = "complex", *args):

        self.wall = wall
        self.pos = pos  # position of image source
        self.sourcePos = sourcePos
        self.sourceStrength = sourceStrength
        self.strength = 1.0
        self.order = order
        self.theta = []  # angle with receiver
        self.r = []  # distance from reciver
        self.img_type = img_type

        if len(args) == 0:
            self.finite_walls = False
        else:
            self.finite_walls = args[0]

    def getRelativeReceiverParameters(self, receiverPos):
        """
        Calculates distance and angle of image source from receiver

        Parameters
        ----------
        receiverPos : microphone position
               3D array with cartesian coordinates

        Returns
        -------
        None.

        """

        for i in range(len(receiverPos)):
            self.r.append(self.pos.getDistance(receiverPos[i]))
            v1 = -self.pos.subtract(self.sourcePos)
            v2 = -self.pos.subtract(receiverPos[i])
            self.theta.append(self.getAngle(v1, v2))

    def calculateAngleLimitsWithOtherWalls(self, wall):
        """
        Calculates limiting angles with the other walls of the room
        for it acts as a virtual source

        Returns
        -------
        angleLims - 2 element array with lower and upper limits

        """

        # vector from source to image source
        v0 = -self.sourcePos.subtract(self.pos)

        # wall edges next to each other should be used, since walls are set
        # clockwise, this is the correct way to do it
        v1 = -self.sourcePos.subtract(wall.posA)
        v2 = -self.sourcePos.subtract(wall.posD)
        angleLims = np.array([-self.getAngle(v0, v1), self.getAngle(v0, v2)])
        return angleLims

    def calculateImageStrength(self, wall, wav_nums):
        """
        Calculates image source strength based on spherical scattering from
        wall

        Parameters
        ----------
        wall     : the wall for which it acts as a virtual source
        wav_nums : array
            wave numbers to evaluate on.

        Returns
        -------
        None

        """

        assert len(wav_nums) == len(wall.wallImpedance)

        Nx = len(wav_nums)
        Ny = len(self.theta)
        [k, r] = np.meshgrid(wav_nums, self.r)
        kr = np.multiply(k, r)

        beta = wall.wallImpedance
        [beta, theta] = np.meshgrid(beta, self.theta)
        gamma0 = np.cos(theta)
        R0 = np.divide(gamma0 - beta, gamma0 + beta)


        if self.img_type == 'complex':
            if self.finite_walls:
                rho = np.divide(gamma0 + beta, np.sqrt(2 * (1 + np.multiply(gamma0, beta))))
                angleLims = self.calculateAngleLimitsWithOtherWalls(wall)
                etaMax = theta - (np.ones((Ny, Nx)) * angleLims[0])
                tMax = -1j * (1 - np.cos(etaMax))

                etaMin = theta - (np.ones((Ny, Nx)) * angleLims[1])
                tMin = -1j * (1 - np.cos(etaMin))

                integral = self.evaluateIntegral(kr, rho, tMax) - self.evaluateIntegral(
                    kr, rho, tMin
                )
                self.strength = R0 + (1 - R0) * (
                    1 - (rho * np.exp(-1j * kr * rho) * integral)
                )

            else:

                R0 = np.divide(gamma0 - beta, gamma0 + beta)
                rho = (
                    1j
                    * kr
                    / 2.0
                    * np.divide(np.power(beta + gamma0, 2), (1 + beta * gamma0))
                )
                F = 1 + 1j * np.sqrt(np.pi * rho) * np.exp(-rho) * erfc(-1j * np.sqrt(rho))
                self.strength = R0 + (1 - R0) * F
        else:
            self.strength = R0

        # attenuate by strength of source
        self.strength *= self.sourceStrength

    def evaluateIntegral(self, c, B, t):

        gammaVal = myGammainc(c, B, t)
        return np.exp(-1j * B * c) * np.divide(
            np.sqrt(c * (t - 1j * B) * gammaVal), np.sqrt(B + 1j * t)
        )

    def getAngle(self, v1, v2):
        """
        Calculates angle between two vectors

        Parameters
        ----------
        v1 : direction vector (3x1)
        v2 : direction vector (3x1)

        Returns
        -------
        Angle in radians

        """
        return np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
