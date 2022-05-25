#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 14:31:41 2021

@author: od0014
"""
import math
import numpy as np
from sympy.utilities.iterables import multiset_permutations


class Point:
    """
    Class that defines a point in 3D cartesian coordinate system
    """
    
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        
    def getDistance(self, p):
        
        """ 
        Returns the Euclidean distance between two 3D positions in 
        cartesian coordinates
        """
        dis = math.sqrt((self.x - p.x) ** 2 + (self.y - p.y) ** 2 + (self.z - p.z) ** 2)
        return dis
    
    
    def subtract(self, p):
        
        return np.array([self.x - p.x, self.y - p.y, self.z - p.z], dtype = float)
    
    def add(self, p):
        
        return Point(self.x + p.x, self.y + p.y, self.z + p.z)
    
    def divide(self, scalar):
        
        return Point(self.x/scalar, self.y/scalar, self.z/scalar)
    
    def equals(self, p, tolerance=0):
        if tolerance == 0:
            if (self.x == p.x) and\
                    (self.y == p.y) and\
                    (self.z == p.z):
                return True
            else:
                return False
        else:
            if self.getDistance(p) < tolerance:
                return True
            else:
                return False

    def getNorm(self):
        
        norm = np.sqrt(self.x * self.x + self.y * self.y + self.z*self.z)
        return norm

    def getDotProduct(self, p):
        
        dotProduct = self.x * p.x + self.y * p.y + self.z * p.z
        return dotProduct

    def avoid_corners(self, wall, safety=1e-2):
        max_x = np.round(np.max([wall.posA.x, wall.posB.x, wall.posC.x, wall.posD.x]), 5)
        min_x = np.round(np.min([wall.posA.x, wall.posB.x, wall.posC.x, wall.posD.x]), 5)

        max_y = np.round(np.max([wall.posA.y, wall.posB.y, wall.posC.y, wall.posD.y]), 5)
        min_y = np.round(np.min([wall.posA.y, wall.posB.y, wall.posC.y, wall.posD.y]), 5)

        max_z = np.round(np.max([wall.posA.z, wall.posB.z, wall.posC.z, wall.posD.z]), 5)
        min_z = np.round(np.min([wall.posA.z, wall.posB.z, wall.posC.z, wall.posD.z]), 5)

        if self.x == min_x:
            s_x = min_x + safety
        elif self.x == max_x:
            s_x = max_x - safety
        else:
            s_x = self.x

        if self.y == min_y:
            s_y = min_y + safety
        elif self.y == max_y:
            s_y = max_y - safety
        else:
            s_y = self.y

        if self.z == min_z:
            s_z = min_z + safety
        elif self.z == max_z:
            s_z = max_z - safety
        else:
            s_z = self.z

        if min_x == max_x:
            s_x = min_x
        if min_y == max_y:
            s_y = min_y
        if min_z == max_z:
            s_z = min_z

        return Point(s_x, s_y, s_z)


class Plane:
    
    """
    Class and helper functions defining a plane in 3D
    """
    
    def __init__(self):
        #plane represented by ax + by + cz + d = 0 and its normal vector

        self.a = 0
        self.b = 0
        self.c = 0
        self.d = 0
        self.normal = np.zeros(3, dtype=float)

    def setPlaneFromPoints(self, posA, posB, posC):
        #posA, posB and posC are 3 points on a plane
        
        #find vector normal to the plane
        arr1 = posB.subtract(posA)
        arr2 = posC.subtract(posA)
        self.normal = np.cross(arr1, arr2)
        
        assert np.dot(self.normal, arr1) == 0.0, "normal vector not right"
        
        self.a = self.normal[0]
        self.b = self.normal[1]
        self.c = self.normal[2]
        
        #scalar component
        self.d = np.dot(-self.normal, [posA.x, posA.y, posA.z])
        # self.d = -(self.a*posB.x + self.b*posB.y + self.x*posB.z)

    def equals(self, otherPlane):
        """check if two planes are equal"""
        
        eps = np.finfo(np.float32).eps
        P = (self.a+eps)/(otherPlane.a+eps)
        Q = (self.b+eps)/(otherPlane.b+eps)
        R = (self.c+eps)/(otherPlane.c+eps)
        S = (self.d+eps)/(otherPlane.d+eps)
        
        if P == Q and Q == R and R == S :
            print(P,Q,R,S)
            return True
        else:
            return False

    def pointInPlane(self, point):
        """check whether Point object is bounded by plane"""

        tol = 1e-12
        
        if abs(self.a*point.x + self.b*point.y + self.c*point.z + self.d) < tol:
            return True
        else:
            return False

    def getPointReflection(self, point):
        """get reflection of a point along the plane"""
        
        # equation of line from (x1,y1,z1) to where it intersects plane
        # (x - x1) / a = (y - y1) / b = (z - z1) / c = k
        # replace x = ak + x1 etc in ax + by + cz + d = 0 to find k
        
        k = -(self.a * point.x + self.b * point.y + self.c * point.z + 
             self.d) / (self.a ** 2 + self.b ** 2 + self.c ** 2)
        
        # where line from point intersects plane, is the midpoint between (x1,y1,z1)
        # and its reflection
        
        refPos = Point(0.0, 0.0, 0.0)
        refPos.x = 2 * (self.a * k + point.x) - point.x
        refPos.y = 2 * (self.b * k + point.y) - point.y
        refPos.z = 2 * (self.c * k + point.z) - point.z
        
        return refPos

    def findLineIntersection(self, posA, posB):
        """find point where a line intersects the plane"""
        
        #two points are enough to define a line
        #equation of a line is (x-x1)/l = (y-y1)/m = (z-z1)/n = k
        
        l = posB.x - posA.x 
        m = posB.y - posA.y
        n = posB.z - posA.z
        
        # replace x with kl + x1 etc and plug into ax + by + cz + d = 0 to find k
        k = -(self.a * posA.x + self.b * posA.y + self.c * posA.z + 
              self.d) / (self.a * l + self.b * m + self.c * n)
        
        # plug in value of k into x = kl+x1 etc to find point of intersection
        
        interPos = Point(0.0, 0.0, 0.0)
        interPos.x = k*l + posA.x
        interPos.y = k*m + posA.y
        interPos.z = k*n + posA.z
        
        assert self.pointInPlane(interPos), "intersection point does not lie on plane!"
        return interPos


class Room:
    """
    Class defining a room with some properties that can be controlled
    """
    def __init__(self):
        self.shape = ''
        self.wallImpedance = list()   #this is a dictionary


class Cuboid(Room):
    """
    Class defining a cuboid room with dimensions and wall positions
    """
    
    def __init__(self, x, y, z):
        self.name = 'cuboid'
        self.nWalls = 6
        self.walls = dict()
        self.wall_tree = WallTree(None,  # parent (WallTree class)
                                  None,  # wall data (Wall class)
                                  None,  # wall type ('floor', ...)
                                  0)  # order
       
        self.x = x
        self.y = y
        self.z = z

    def setWallPosition(self, refOrder, wallImpedances):
        ind = [-1, 1]
        
        # loop through reflection orders and create virtual reflection walls
        for i in range(1, refOrder+1):
            # first order reflectors are the actual room walls
            if i == 1:
                #create empty lists - this is going to be a dictionary of lists
                self.walls['floor'] = list()
                self.walls['ceiling'] = list()
                self.walls['left'] = list()
                self.walls['right'] = list()
                self.walls['front'] = list()
                self.walls['back'] = list()

                # append to empty list - these are added clockwise for a reason
                self.walls['floor'].append(Wall(Point(0,0,0), Point(self.x, 0, 0),
                                                Point(self.x, 0, self.z), Point(0, 0, self.z), i))
                
                self.walls['ceiling'].append(Wall(Point(0, self.y, 0), Point(0, self.y, self.z),
                                                  Point(self.x, self.y, self.z), Point(self.x, self.y, 0), i))
                
                self.walls['left'].append(Wall(Point(0,0,0), Point(0, self.y, 0), 
                                               Point(0, self.y, self.z), Point(0, 0, self.z), i))
                
                self.walls['right'].append(Wall(Point(self.x, 0, 0), Point(self.x, self.y, 0),
                                                Point(self.x, self.y, self.z), Point(self.x, 0, self.z), i))
                
                self.walls['front'].append(Wall(Point(0,0,self.z), Point(0, self.y, self.z), 
                                                Point(self.x, self.y, self.z), Point(self.x, 0, self.z), i))
                
                self.walls['back'].append(Wall(Point(0,0,0), Point(0, self.y, 0), 
                                               Point(self.x, self.y, 0), Point(self.x, 0, 0), i))

            else:
                # iterate through the 6 walls
                for key in self.walls.keys():
                    # for ith order reflection - for each wall (N-1)^{i-1} extra wall reflectors should be added.
                    # reflectedWalls should be the walls added in the last reflection order
                    reflectedWalls = [wall for wall in self.walls[key] if wall.order == (i-1)]
                    
                    for originalWall in reflectedWalls:
                        # k = [-1, 1]
                        for k in ind:
                            # translate x, y and z from original wall positions
                            
                            wallX = Wall(originalWall.posA.add(Point(k*self.x, 0, 0)),
                                         originalWall.posB.add(Point(k*self.x, 0, 0)),
                                         originalWall.posC.add(Point(k*self.x, 0, 0)),
                                         originalWall.posD.add(Point(k*self.x, 0, 0)),
                                         i)

                            wallY = Wall(originalWall.posA.add(Point(0, k*self.y, 0)),
                                         originalWall.posB.add(Point(0, k*self.y, 0)),
                                         originalWall.posC.add(Point(0, k*self.y, 0)),
                                         originalWall.posD.add(Point(0, k*self.y, 0)),
                                         i)

                            wallZ = Wall(originalWall.posA.add(Point(0, 0, k*self.z)),
                                         originalWall.posB.add(Point(0, 0, k*self.z)),
                                         originalWall.posC.add(Point(0, 0, k*self.z)),
                                         originalWall.posD.add(Point(0, 0, k*self.z)),
                                         i)

                            # Create a list of all values in list of dictionaries -
                            list_of_all_walls = [wall for key in self.walls for wall in self.walls[key]] 
                                                 # if wall.order == i]
                            # print("Number of reflecting walls is " + str(len(list_of_all_walls)))

                            # if reflector already exists in dictionary, dont add it
                            xFlag = True
                            yFlag = True
                            zFlag = True
                            
                            for wall in list_of_all_walls:
                                if wallX.equals(wall):
                                    xFlag = False
                                if wallY.equals(wall):
                                    yFlag = False
                                if wallZ.equals(wall):
                                    zFlag = False

                            if xFlag:
                                self.walls[key].append(wallX)
                            if yFlag:
                                self.walls[key].append(wallY)
                            if zFlag:
                                self.walls[key].append(wallZ)

        # now that we have added all the walls, set the plane coefficients and the wall filters
        k = 0
        for key in self.walls.keys():
            #add wall filters
            self.walls[key][0].setWallImpedance(wallImpedances[k])
            k += 1
            
            # set plane coefficients
            for item in self.walls[key]:
                item.setPlaneCoefficients()

        return self.walls



    def setWallPositionWithWallTree(self, source_position, refOrder, wallImpedances):
        self.wall_tree.img_source_position = source_position

        # loop through reflection orders and create virtual reflection walls
        for i in range(1, refOrder+1):
            # first order reflectors are the actual room walls
            if i == 1:
                # create structs with wall data, children lists will be filled later
                self.wall_tree.children.append(WallTree(self.wall_tree,  # parent
                                                        Wall(Point(0, 0, 0),
                                                             Point(self.x, 0, 0),
                                                             Point(self.x, 0, self.z),
                                                             Point(0, 0, self.z),
                                                             i),  # wall data
                                                        'floor',  # wall type
                                                        i))  # order
                self.wall_tree.children.append(WallTree(self.wall_tree,  # parent
                                                        Wall(Point(0, self.y, 0),
                                                             Point(0, self.y, self.z),
                                                             Point(self.x, self.y, self.z),
                                                             Point(self.x, self.y, 0),
                                                             i),  # wall data
                                                        'ceiling',  # wall type
                                                        i))  # order
                self.wall_tree.children.append(WallTree(self.wall_tree,  # parent
                                                        Wall(Point(0, 0, 0),
                                                             Point(0, self.y, 0),
                                                             Point(0, self.y, self.z),
                                                             Point(0, 0, self.z),
                                                             i),  # wall data
                                                        'left',  # wall type
                                                        i))  # order
                self.wall_tree.children.append(WallTree(self.wall_tree,  # parent
                                                        Wall(Point(self.x, 0, 0),
                                                             Point(self.x, self.y, 0),
                                                             Point(self.x, self.y, self.z),
                                                             Point(self.x, 0, self.z),
                                                             i),  # wall data
                                                        'right',  # wall type
                                                        i))  # order
                self.wall_tree.children.append(WallTree(self.wall_tree,  # parent
                                                        Wall(Point(0, 0, self.z),
                                                             Point(0, self.y, self.z),
                                                             Point(self.x, self.y, self.z),
                                                             Point(self.x, 0, self.z),
                                                             i),  # wall data
                                                        'front',  # wall type
                                                        i))  # order
                self.wall_tree.children.append(WallTree(self.wall_tree,  # parent
                                                        Wall(Point(0, 0, 0),
                                                             Point(0, self.y, 0),
                                                             Point(self.x, self.y, 0),
                                                             Point(self.x, 0, 0),
                                                             i),  # wall data
                                                        'back',  # wall type
                                                        i))  # order

                for k, wall in enumerate(self.wall_tree.children):
                    # add wall filters
                    wall.wall_data.setWallImpedance(wallImpedances[k])
                    # set plane coefficients
                    wall.wall_data.setPlaneCoefficients()

                    # compute image source locations
                    wall.compute_img_source()
            else:
                # for ith order reflection - for each wall (N-1)^{i-1} extra wall reflectors should be added.
                # reflectedWalls are the walls added in the previous reflection order
                reflected_walls = self.wall_tree.get_children_of_order(i-1)

                for reference_wall in reflected_walls:
                    reference_plane = reference_wall.wall_data.plane
                    # take all the walls in the same "imaginary room" as this wall,
                    # excluding this wall itself
                    siblings = [wall for wall in reference_wall.parent.children if wall != reference_wall]
                    # but including this wall's parent, if it is a valid wall
                    if reference_wall.parent.wall_data is not None:
                        siblings = siblings + [reference_wall.parent]

                    # each of these walls is mirrored w.r.t. the current wall to create a mirrored room
                    for sibling in siblings:
                        newA = reference_plane.getPointReflection(sibling.wall_data.posA)
                        newB = reference_plane.getPointReflection(sibling.wall_data.posB)
                        newC = reference_plane.getPointReflection(sibling.wall_data.posC)
                        newD = reference_plane.getPointReflection(sibling.wall_data.posD)

                        new_reflection = Wall(newA, newB, newC, newD, i)
                        # this should have the impedance of wall with same wall_type of previous order
                        new_reflection.setWallImpedance(sibling.wall_data.wallImpedance)   
                        new_reflection.setPlaneCoefficients()
                        
# =============================================================================
#                         Because of this mirroring, a wall that is ceiling is 
#                         incorrectly labeled as floor
#                         In general, opposite walls are incorrectly labeled. 
#                         So, to fix that
# =============================================================================
                        if sibling.wall_type == self.getOppositeWallLabel(reference_wall.wall_type):
                            new_leaf = new_leaf = WallTree(reference_wall,  # parent
                                            new_reflection,  # wall data
                                            self.getOppositeWallLabel(sibling.wall_type),  # wall type
                                            i)  # order
                        else:
                            new_leaf = WallTree(reference_wall,  # parent
                                            new_reflection,  # wall data
                                            sibling.wall_type,  # wall type
                                            i)  # order
                        # compute image source location
                        new_leaf.compute_img_source()

                        # TODO: check visibility of entire wall before appending to tree?
                        reference_wall.children.append(new_leaf)

        return self.wall_tree
    
    
    def getOppositeWallLabel(self,wall_type):

        if wall_type == 'floor':
            return 'ceiling'
        elif wall_type == 'ceiling':
            return  'floor'
        elif wall_type == 'left':
            return 'right'
        elif wall_type == 'right':
            return 'left'
        elif wall_type == 'front':
            return 'back'
        elif wall_type == 'back':
            return 'front'



class Wall:
    
    """
    Class defining a wall in a room, which represents a plane in 3D 
    """
    
    def __init__(self, posA, posB, posC, posD, order):
        self.posA = posA
        self.posB = posB
        self.posC = posC
        self.posD = posD
        self.plane = Plane()
        self.wallImpedance = list()
        self.order = order

    def setWallImpedance(self, impedance):
        #list of impedances at desired wave numbers
        self.wallImpedance = impedance


    def setPlaneCoefficients(self):
        self.plane.setPlaneFromPoints(self.posA, self.posB, self.posC)

    def equals(self, otherWall):
        """checks if two walls are identical"""
     
        A = np.array([np.array([self.posA.x, self.posA.y, self.posA.z]), \
                      np.array([self.posB.x, self.posB.y, self.posB.z]), \
                      np.array([self.posC.x, self.posC.y, self.posC.z]), \
                          np.array([self.posD.x, self.posD.y, self.posD.z])])

        B = np.array([np.array([otherWall.posA.x, otherWall.posA.y, otherWall.posA.z]),\
                     np.array([otherWall.posB.x, otherWall.posB.y, otherWall.posB.z]),\
                     np.array([otherWall.posC.x, otherWall.posC.y, otherWall.posC.z]),\
                     np.array([otherWall.posD.x, otherWall.posD.y, otherWall.posD.z])])                                           
        
            
        # any permutation of these points (matrix rows) should be equal for the walls to be same
        A_perm = np.empty(np.shape(A)) 
        row_pos = [0,1,2,3]
        permute = multiset_permutations(row_pos)
       
        for p in permute:
            A_perm[:,:] = A[p,:]
          
            if np.array_equal(A_perm,B):

                return True
            
        return False

    def boundedByRoom(self, point, wallType, room):
        """checks whether an image source is bounded by the room the wall is located in"""
        
        midPoint = self.posA.add(self.posC)
        midPoint = midPoint.divide(2)

        # project vector connecting plane mid-point and external point and check
        # if it is less than the corresponding room dimension
        
        vec = midPoint.subtract(point)
        unit_vec = vec/ np.linalg.norm(vec)
        
        normal = np.abs(self.plane.normal)
        unit_normal = normal/np.linalg.norm(normal)
        
        angle = np.arccos(np.dot(unit_vec, unit_normal))
        projection = vec * np.cos(angle)

        projection_x = projection[0]
        projection_y = projection[1]
        projection_z = projection[2]


        if wallType == 'floor' or wallType == 'ceiling':
            if point.x > (midPoint.x - room.shape.x/2) and point.x < (midPoint.x + room.shape.x/2) \
                and point.z > (midPoint.z - room.shape.z/2) and (point.z < midPoint.z + room.shape.z/2) \
                    and abs(projection_y) < room.shape.y:
                        # print('Wall type ' + wallType)
                        # print('Wall midpoint =', midPoint.x, midPoint.y, midPoint.z)
                        # print('Image location =', point.x, point.y, point.z)
                        return True
            else:
                return False
        
        elif wallType == 'left' or wallType == 'right':
            if point.y > midPoint.y - room.shape.y/2 and point.y < midPoint.y+ room.shape.y/2 \
                and point.z > midPoint.z - room.shape.z/2 and point.z < midPoint.z + room.shape.z/2 \
                    and  abs(projection_x) < room.shape.x:
                        return True
            else:
                return False

        elif wallType == 'front' or wallType == 'back':
            if point.x > midPoint.x - room.shape.x/2 and point.x < midPoint.x + room.shape.x/2 \
                and point.y > midPoint.y - room.shape.y/2 and point.y < midPoint.y + room.shape.y/2 \
                    and abs(projection_z) < room.shape.z:
                        return True
            else:
                return False

   

    def boundedByWall(self, p):
        """checks if a point is bounded by the 4 points defining the wall"""

        p_x = np.round(p.x, 5)
        p_y = np.round(p.y, 5)
        p_z = np.round(p.z, 5)
        
        max_x = np.round(np.max([self.posA.x, self.posB.x, self.posC.x, self.posD.x]), 5)
        min_x = np.round(np.min([self.posA.x, self.posB.x, self.posC.x, self.posD.x]), 5)
        
        max_y = np.round(np.max([self.posA.y, self.posB.y, self.posC.y, self.posD.y]), 5)
        min_y = np.round(np.min([self.posA.y, self.posB.y, self.posC.y, self.posD.y]), 5)
        
        max_z = np.round(np.max([self.posA.z, self.posB.z, self.posC.z, self.posD.z]), 5)
        min_z = np.round(np.min([self.posA.z, self.posB.z, self.posC.z, self.posD.z]), 5)
        
        if (p_x < min_x or p_x > max_x or
            p_y < min_y or p_y > max_y or
            p_z < min_z or p_z > max_z):
            return False
        else:
            return True




class WallTree:

    """
    Class that defines the imaginary wall structure,
    linking each imaginary wall to the parent (last wall this is "seen through")
    and children (first walls that are "seen through" this
    """

    def __init__(self, parent, wall_data, wall_type, order):
        self.parent = parent
        self.wall_data = wall_data
        self.wall_type = wall_type
        self.order = order
        self.children = list()
        self.img_source_position = None

    def get_real_counterpart(self):
        root = self
        while root.order > 0:
            root = root.parent

        real_walls = root.children

        for real_wall in real_walls:
            if real_wall.wall_type == self.wall_type:
                return real_wall

        raise Exception("Could not find the real counterpart to a wall")

    def compute_img_source(self):
        if self.parent.img_source_position is None:
            self.parent.compute_img_source()
        self.img_source_position = self.wall_data.plane.getPointReflection(self.parent.img_source_position)

        return self.img_source_position

    def get_all_children(self):
        recursive_list = sum([child.get_all_children() for child in self.children], [])
        return [self] + recursive_list

    def get_children_of_order(self, order):
        # note: if no nodes of the correct order exist, an empty list will bre returned

        if order == self.order:
            # no need to look further into children, they have higher orders
            return [self]

        elif order > self.order:
            recursive_list = sum([child.get_children_of_order(order) for child in self.children], [])
            return recursive_list

        else:  # order < self.order
            raise("Either there's an incorrect 'order' value, " +
                  "or the function was called at the wrong node.")
