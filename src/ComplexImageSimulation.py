#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 11:36:46 2022

@author: od0014
"""

import numpy as np
import ComplexImageMethod.src.Reflection as ref



class ComplexImageSimulation:
    
    def __init__(self, room, source, mic_array, wave_num, order, *args):
        self.room = room
        self.source = source
        self.mic_array = mic_array
        self.wave_num = wave_num
        self.maxOrder = order      #reflection order
        if len(args) == 0:
            self.finiteWall = False
        else:
            self.finiteWall = args[0]
        
                

        
    def run(self):
               
        #list of real and virtual walls upto given order
        wall_list = self.room.shape.setWallPositionWithWallTree(self.source, self.maxOrder, self.room.wallImpedance)

        reflectedPressure = np.zeros((len(self.mic_array), len(self.wave_num)), dtype = complex)
        
        #create reflection object
        reflect = ref.Reflection(self.source, self.maxOrder)
        
        real_walls = wall_list.children
        
        #loop through reflection orders
        for order in range(1, self.maxOrder+1):
            
            # loop through each plane - floor, ceiling etc
            for mainWall in real_walls:
            
                #walls corresponding to this reflection order
                whichWall = mainWall.wall_type
                walls_thisOrder = [wall for wall in wall_list.get_children_of_order(order) if wall.wall_type == whichWall]
                
                
                #find image source corresponding to particular order and linked with thisWall
                for thisWall in walls_thisOrder:

                    curImageSource = reflect.createImageSource(thisWall, whichWall, order, self.room, self.finiteWall)
                    
                    #find receiver distance and angle from image source
                    curImageSource.getRelativeReceiverParameters(self.mic_array)
                    
                    curImageSource.calculateImageStrength(thisWall.wall_data, whichWall, self.wave_num)
                    
                    #this acts as virtual sources for other walls
                    for otherWall in real_walls: 
                        
                        wallName = otherWall.wall_type
                        
                        curImageSource.calculateImageStrength(otherWall.wall_data, wallName, self.wave_num)

                        if wallName == whichWall:
                             
                            # calculate reflected pressure due to this image source
                            [k, r] = np.meshgrid(self.wave_num, curImageSource.r)
                            kr = np.multiply(k,r)                     
                            reflectedPressure += np.divide(np.exp(1j*kr),kr) * curImageSource.strength[wallName]

                    
                    # now is the right time to append it to the list of image sources
                    reflect.addImageSource(curImageSource, order)
                            
                            
            print('Reflected pressure calculated for reflection order ', str(order))
           
        
        
        #calculate direct pressure
        dis = np.zeros(len(self.mic_array))
        for i in range(len(self.mic_array)):
            dis[i] = (self.source.getDistance(self.mic_array[i]))
        
        [k, r0] = np.meshgrid(self.wave_num, dis)
        kr0 = np.multiply(k,r0) 
        directPressure = np.divide(np.exp(1j*kr0),kr0) 
        
        #total pressure is the sum of direct and reflceted pressures
        totalPressure = directPressure + reflectedPressure
        
        
        return reflectedPressure, totalPressure
            
                            
                    
                    

        
        
        
        
        


        
        
        
    
       
