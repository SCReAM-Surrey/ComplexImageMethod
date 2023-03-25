#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 18:22:27 2022

@author: od0014
"""

import math
import numpy as np


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

        return np.array([self.x - p.x, self.y - p.y, self.z - p.z], dtype=float)

    def add(self, p):

        return Point(self.x + p.x, self.y + p.y, self.z + p.z)

    def divide(self, scalar):

        return Point(self.x / scalar, self.y / scalar, self.z / scalar)

    def equals(self, p, tolerance=0):
        if tolerance == 0:
            if (self.x == p.x) and (self.y == p.y) and (self.z == p.z):
                return True
            else:
                return False
        else:
            if self.getDistance(p) < tolerance:
                return True
            else:
                return False

    def getNorm(self):

        norm = np.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)
        return norm

    def getDotProduct(self, p):

        dotProduct = self.x * p.x + self.y * p.y + self.z * p.z
        return dotProduct


class Plane:

    """
    Class and helper functions defining a plane in 3D
    """

    def __init__(self, a, b, c, d):
        # plane represented by ax + by + cz + d = 0 and its normal vector

        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.normal = np.array([a, b, c], dtype=float)

    def setPlaneFromPoints(self, posA, posB, posC):
        # posA, posB and posC are 3 points on a plane

        # find vector normal to the plane
        arr1 = posB.subtract(posA)
        arr2 = posC.subtract(posA)
        self.normal = np.cross(arr1, arr2)

        assert np.dot(self.normal, arr1) == 0.0, "normal vector not right"

        self.a = self.normal[0]
        self.b = self.normal[1]
        self.c = self.normal[2]

        # scalar component
        self.d = np.dot(-self.normal, [posA.x, posA.y, posA.z])
        # self.d = -(self.a*posB.x + self.b*posB.y + self.x*posB.z)

    def equals(self, otherPlane):
        """check if two planes are equal"""

        eps = np.finfo(np.float32).eps
        P = (self.a + eps) / (otherPlane.a + eps)
        Q = (self.b + eps) / (otherPlane.b + eps)
        R = (self.c + eps) / (otherPlane.c + eps)
        S = (self.d + eps) / (otherPlane.d + eps)

        if P == Q and Q == R and R == S:
            print(P, Q, R, S)
            return True
        else:
            return False

    def pointInPlane(self, point):
        """check whether Point object is bounded by plane"""

        tol = 1e-12

        if abs(self.a * point.x + self.b * point.y + self.c * point.z + self.d) < tol:
            return True
        else:
            return False

    def getPointReflection(self, point):
        """get reflection of a point along the plane"""

        # equation of line from (x1,y1,z1) to where it intersects plane
        # (x - x1) / a = (y - y1) / b = (z - z1) / c = k
        # replace x = ak + x1 etc in ax + by + cz + d = 0 to find k

        k = -(self.a * point.x + self.b * point.y + self.c * point.z + self.d) / (
            self.a**2 + self.b**2 + self.c**2
        )

        # where line from point intersects plane, is the midpoint between (x1,y1,z1)
        # and its reflection

        refPos = Point(0.0, 0.0, 0.0)
        refPos.x = 2 * (self.a * k + point.x) - point.x
        refPos.y = 2 * (self.b * k + point.y) - point.y
        refPos.z = 2 * (self.c * k + point.z) - point.z

        return refPos

    def findLineIntersection(self, posA, posB):
        """find point where a line intersects the plane"""

        # two points are enough to define a line
        # equation of a line is (x-x1)/l = (y-y1)/m = (z-z1)/n = k

        l = posB.x - posA.x
        m = posB.y - posA.y
        n = posB.z - posA.z

        # replace x with kl + x1 etc and plug into ax + by + cz + d = 0 to find k
        k = -(self.a * posA.x + self.b * posA.y + self.c * posA.z + self.d) / (
            self.a * l + self.b * m + self.c * n
        )

        # plug in value of k into x = kl+x1 etc to find point of intersection

        interPos = Point(0.0, 0.0, 0.0)
        interPos.x = k * l + posA.x
        interPos.y = k * m + posA.y
        interPos.z = k * n + posA.z

        assert self.pointInPlane(interPos), "intersection point does not lie on plane!"
        return interPos


class Room:
    """
    Class defining a room with some properties that can be controlled
    """

    def __init__(self):
        self.shape = ""
        self.wallImpedance = dict()  # this is a dictionary


class Cuboid(Room):
    """
    Class defining a cuboid room with dimensions and wall positions
    """

    def __init__(self, x, y, z):
        self.name = "cuboid"
        self.nWalls = 6
        self.walls = dict()

        self.x = x
        self.y = y
        self.z = z

    def setWallPosition(self, wallImpedances):

        self.walls["floor"] = Wall(
            "floor", self.x, self.y, self.z, wallImpedances["floor"]
        )
        self.walls["ceiling"] = Wall(
            "ceiling", self.x, self.y, self.z, wallImpedances["ceiling"]
        )
        self.walls["left"] = Wall(
            "left", self.x, self.y, self.z, wallImpedances["left"]
        )
        self.walls["right"] = Wall(
            "right", self.x, self.y, self.z, wallImpedances["right"]
        )
        self.walls["front"] = Wall(
            "front", self.x, self.y, self.z, wallImpedances["front"]
        )
        self.walls["back"] = Wall(
            "back", self.x, self.y, self.z, wallImpedances["back"]
        )

        return self.walls


class Wall:

    """
    Class defining a wall in a room, which represents a plane in 3D
    z   y
    |  /
    | /
    |/
    ----------- x
    """

    def __init__(self, wallID, L, W, H, impedance):

        self.name = wallID
        # these are infinite planes,
        # I am setting it such that the normal points outwards
        if self.name == "floor":
            self.plane = Plane(0, 0, -1, 0)
            self.posA = Point(0, 0, 0)
            self.posD = Point(L, 0, 0)
        elif self.name == "ceiling":
            self.plane = Plane(0, 0, 1, -H)
            self.posA = Point(0, 0, H)
            self.posD = Point(L, 0, H)
        elif self.name == "left":
            self.plane = Plane(-1, 0, 0, 0)
            self.posA = Point(0, 0, 0)
            self.posD = Point(0, W, 0)
        elif self.name == "right":
            self.plane = Plane(1, 0, 0, -L)
            self.posA = Point(L, 0, 0)
            self.posD = Point(L, W, 0)
        elif self.name == "front":
            self.plane = Plane(0, -1, 0, 0)
            self.posA = Point(0, 0, 0)
            self.posD = Point(L, 0, 0)
        elif self.name == "back":
            self.plane = Plane(0, 1, 0, -W)
            self.posA = Point(0, W, 0)
            self.posD = Point(L, W, 0)

        # list of impedances at desired wave numbers
        self.wallImpedance = impedance

    def isPointBehindWall(self, point):
        distance = np.dot(self.plane.normal, self.posA.subtract(point))
        if distance < 0:
            return True
        else:
            return False
