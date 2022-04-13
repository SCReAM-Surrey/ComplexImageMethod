#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 11:27:59 2022

@author: od0014
"""

import ComplexImageMethod.src.ImageSource as ImgSrc



class Reflection:
    """
    Class that determines wall node positions by calculating 1st order reflections
    """
    
    def __init__(self, sourcePos):
        self.sourcePos = sourcePos
        self.imageSrc = list()
        self.roomWalls = dict()     # keeps track of the physical walls in the room

    
    def getReflectionAlongPlane(self, objPos, wall):
        return wall.plane.getPointReflection(objPos)

    
    def createImageSource(self, wall, wallType, order, room, finiteFlag):
        # find the correct image source to reflect along the wall

        if order == 1:
            self.roomWalls[wallType] = wall
            srcPos = self.sourcePos
            srcStrength = {
                "floor": 1.0,
                "ceiling": 1.0,
                "left": 1.0,
                "right": 1.0,
                "front": 1.0,
                "back": 1.0}
        else:
            for img in self.imageSrc:
                # look only at image sources of previous order
                if img.order == order - 1:
                    # check if image-source is bounded by the wall
                    if wall.boundedByRoom(img.pos, wallType, room):
                        srcPos = img.pos
                        srcStrength = img.strength
                        break
                else:
                    continue


        # add to existing image sources
        curImageSource = ImgSrc.ImageSource(srcPos, srcStrength, wall, order, finiteFlag)
        return curImageSource
        
        
    def addImageSource(self, imgSource):
        self.imageSrc.append(imgSource)


                

        


   