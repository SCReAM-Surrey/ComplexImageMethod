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
    
    def __init__(self, sourcePos, maxOrder):
        self.sourcePos = sourcePos
        self.imageSrc = {k: [] for k in range(1, maxOrder+1)}  # dictionary of lists, sorted by order
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
            # look only at image sources of previous order
            for img in self.imageSrc[order-1]:
                if wall.wall_data.boundedByRoom(img.pos, wallType, room):
                    srcPos = img.pos
                    srcStrength = img.strength
                    # if wallType == img.wall.wall_type:
                    #     oppImg = self.getImageOnOppositeWall(img,order)
                    #     srcStrength = oppImg.strength
                    # else:
                    #     srcStrength = img.strength
                    break

        curImageSource = ImgSrc.ImageSource(srcPos, srcStrength, wall, order, finiteFlag)
        return curImageSource
        
        
    def addImageSource(self, imgSource, order):
        # add to existing image sources
        self.imageSrc[order].append(imgSource)
        
    
            
    def getImageOnOppositeWall(self,img,order):
        # I don't know how to explain this but it seemed necessary to me
        # to get the right source
        if img.wall.wall_type == 'floor':
            oppositeWall = 'ceiling'
        elif img.wall.wall_type == 'ceiling':
            oppositeWall = 'floor'
        elif img.wall.wall_type == 'left':
            oppositeWall = 'right'
        elif img.wall.wall_type == 'right':
            oppositeWall = 'left'
        elif img.wall.wall_type == 'front':
            oppositeWall = 'back'
        elif img.wall.wall_type == 'back':
            oppositeWall = 'front'
            
        oppositeImage = [image for image in self.imageSrc[order-1] if 
                             image.wall.wall_type == oppositeWall]
        return oppositeImage[0]
    
    
   
        


                

        


   