#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 15:12:36 2021

@author: od0014
"""

import numpy as np
import mpmath as mp
from scipy.special import erfc



def myWhittaker(rho):
    ypts = []
    rho_list = rho.tolist()
    for i in range(len(rho_list)):
        for x in rho_list[i]:
            ypts.append(mp.whitw(mp.mpf(-0.25), mp.mpf(0.25), x))
    return np.reshape(np.array(ypts, dtype=complex), np.shape(rho))


def myGammainc(c,B,t):
    ypts = []
    B_list = B.tolist()
    c_list = c.to_list()
    t_list = t.to_list()
    for i in range(len(c_list)):
        for k in range(len(c_list[i])):
            ypts.append(mp.gammainc(0.5, c_list[i][k]*(t_list[i][k]-1j*B_list[i][k])))
            
    return np.reshape(np.array(ypts, dtype=complex), np.shape(c))
            
    


class ImageSource:

    """
    Class that keeps track of the image sources and their orders
    """

    def __init__(self, sourcePos, sourceStrength, wall, order, *args):
        
        self.sourcePos = sourcePos
        self.sourceStrength = sourceStrength
        self.order = order
        self.wall = wall
        # self.pos = wall.plane.getPointReflection(sourcePos)     #position of image source
        self.pos = wall.compute_img_source()         #position of image source
        self.strength = dict()
        self.theta = []      #angle with receiver
        self.r = []         #distance from reciver
        # self.otherWalls = dict()
        # self.angleLims = []   #angle between IS wall edges
        
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
        

        #vector from source to image source        
        v0 = -self.sourcePos.subtract(self.pos)
        
        #wall edges next to each other should be used, since walls are set
        #clockwise, this is the correct way to do it
        v1 = -self.sourcePos.subtract(wall.posA)
        v2 = -self.sourcePos.subtract(wall.posD)
        angleLims = np.array([-self.getAngle(v0,v1), self.getAngle(v0,v2)])
        return angleLims
    
    
    def calculateImageStrength(self, wall, wallID, wav_nums):
        """
        Calculates image source strength based on spherical scattering from 
        wall

        Parameters
        ----------
        wall     : the wall for which it acts as a virtual source
        wallID   : type of wall - floor, ceiling etc
        wav_nums : array
            wave numbers to evaluate on.

        Returns
        -------
        None

        """
        
        assert len(wav_nums) == len(wall.wallImpedance)
        angleLims = self.calculateAngleLimitsWithOtherWalls(wall)
        
        
        Nx = len(wav_nums)
        Ny= len(self.theta)
        [k, r] = np.meshgrid(wav_nums, self.r)
        kr = np.multiply(k,r)


        beta = wall.wallImpedance
        [beta, theta] = np.meshgrid(beta, self.theta)
        gamma0 = np.cos(theta)


        
        if self.finite_walls:
            R0 = np.divide(gamma0 - beta, gamma0 + beta)
            rho = np.divide(gamma0 + beta, np.sqrt(2*(1 + np.multiply(gamma0, beta))))
            etaMax = theta - (np.ones((Ny,Nx)) * angleLims[0])
            tMax = -1j * (1 - np.cos(etaMax))
        
            etaMin = theta - (np.ones((Ny,Nx)) * angleLims[1])
            tMin = -1j * (1-np.cos(etaMin))
            
            integral = self.evaluateIntegral(kr, rho, tMax) - self.evaluateIntegral(kr, rho, tMin) 
            self.strength[wallID] = R0 + (1-R0)*(1 - (rho*np.exp(-1j*kr*rho)*integral))
            
            
        else:
            
            
            # slower version with mpmath
            # rho = -1j*np.multiply(kr, np.power(rho,2))
            # self.strength[wallID] = R0 + (1-R0)*(1 - np.sqrt(np.pi*rho)*np.exp(rho)*(1 - np.sqrt(erfc(rho))))

            # faster version with scipy
            R0 = np.divide(gamma0*beta - 1.0, beta*gamma0 + 1.0)
            # R0 = np.divide(gamma0 - beta, gamma0 + beta)

            w = np.sqrt(1j*kr/2.0) * (beta+gamma0) 
            self.strength[wallID] = R0 + (1-R0)*(1 + (1j*w*np.sqrt(np.pi)*np.exp(-np.power(w,2))*
                                                      erfc(-1j*w)))
            
        #attenuate by strength of source
        self.strength[wallID] *= self.sourceStrength[wallID]
      
            
            
        
    def evaluateIntegral(self,c,B,t):
        
        gammaVal = myGammainc(c,B,t)
        return np.exp(-1j*B*c)*np.divide(np.sqrt(c*(t-1j*B)*gammaVal), np.sqrt(B+1j*t))
            
       
        
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
        return  np.arccos(np.dot(v1, v2)/(np.linalg.norm(v1) * np.linalg.norm(v2)))
    
    
        

        
        















