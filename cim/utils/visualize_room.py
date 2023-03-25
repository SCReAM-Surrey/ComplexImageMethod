#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 16:25:19 2021

@author: od0014
Some of this code is to be found online.
Helper class to visualize the image method
"""
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


def cuboid_data(center, size):
    """
    Create a data array for cuboid plotting.


    ============= ================================================
    Argument      Description
    ============= ================================================
    center        center of the cuboid, triple
    size          size of the cuboid, triple, (x_depth,y_width,z_height)
    :type size: tuple, numpy.array, list
    :param size: size of the cuboid, triple, (x_depth,y_width,z_height)
    :type center: tuple, numpy.array, list
    :param center: center of the cuboid, triple, (x,y,z)


    """

    # suppose axis direction: x: to left; y: to inside; z: to upper
    # get the (left, outside, bottom) point
    o = [a - b / 2 for a, b in zip(center, size)]
    # get the depth, width, and height
    l, w, h = size
    x = np.array(
        [
            [
                o[0],
                o[0] + l,
                o[0] + l,
                o[0],
                o[0],
            ],  # x coordinate of points in bottom surface
            [
                o[0],
                o[0] + l,
                o[0] + l,
                o[0],
                o[0],
            ],  # x coordinate of points in upper surface
            [
                o[0],
                o[0] + l,
                o[0] + l,
                o[0],
                o[0],
            ],  # x coordinate of points in outside surface
            [o[0], o[0] + l, o[0] + l, o[0], o[0]],
        ]
    )  # x coordinate of points in inside surface

    y = np.array(
        [
            [
                o[1],
                o[1],
                o[1] + w,
                o[1] + w,
                o[1],
            ],  # y coordinate of points in bottom surface
            [
                o[1],
                o[1],
                o[1] + w,
                o[1] + w,
                o[1],
            ],  # y coordinate of points in upper surface
            [o[1], o[1], o[1], o[1], o[1]],  # y coordinate of points in outside surface
            [o[1] + w, o[1] + w, o[1] + w, o[1] + w, o[1] + w],
        ]
    )  # y coordinate of points in inside surface

    z = np.array(
        [
            [o[2], o[2], o[2], o[2], o[2]],  # z coordinate of points in bottom surface
            [
                o[2] + h,
                o[2] + h,
                o[2] + h,
                o[2] + h,
                o[2] + h,
            ],  # z coordinate of points in upper surface
            [
                o[2],
                o[2],
                o[2] + h,
                o[2] + h,
                o[2],
            ],  # z coordinate of points in outside surface
            [o[2], o[2], o[2] + h, o[2] + h, o[2]],
        ]
    )  # z coordinate of points in inside surface

    return x, y, z


def plot_room(ax, center, depth, width, height):
    X, Y, Z = cuboid_data(center, (depth, width, height))
    ax.plot_surface(X, Y, Z, color="b", rstride=1, cstride=1, alpha=0.1)
    ax.set_xlabel("X")
    ax.set_xlim(center[0] - depth / 2, center[0] + depth / 2)
    ax.set_ylabel("Y")
    ax.set_ylim(center[1] - width / 2, center[1] + width / 2)
    ax.set_zlabel("Z")
    ax.set_zlim(center[2] - height / 2, center[2] + height / 2)


def plot_plane(ax, wall, wallType):

    midpoint = wall.posA.add(wall.posC)
    midpoint = midpoint.divide(2)

    dim = np.abs(wall.posA.subtract(wall.posB))
    depth = int(dim[0])
    width = int(dim[1])
    height = int(dim[2])

    if wallType == "front" or wallType == "back":
        xx, yy = np.meshgrid(range(depth), range(width))

        # calculate corresponding z
        zz = (
            (-wall.plane.normal[0] * xx - wall.plane.normal[1] * yy - wall.plane.d)
            * 1.0
            / wall.plane.normal[2]
        )

    elif wallType == "floor" or wallType == "ceiling":
        xx, zz = np.meshgrid(range(depth), range(height))

        yy = (
            (-wall.plane.normal[0] * xx - wall.plane.normal[1] * zz - wall.plane.d)
            * 1.0
            / wall.plane.normal[2]
        )

    elif wallType == "left" or wallType == "right":
        yy, zz = np.meshgrid(range(width), range(height))

        xx = (
            (-wall.plane.normal[0] * yy - wall.plane.normal[1] * zz - wall.plane.d)
            * 1.0
            / wall.plane.normal[2]
        )

    ax.plot_surface(xx, yy, zz, color="r", rstride=1, cstride=1, alpha=0.2)


def plot_point(ax, point, type, *args):

    if type == "node":
        mark = "o"
        order = args[0] - 1
        col = ["tab:blue", "tab:red", "tab:orange", "tab:purple", "tab:cyan"]
    elif type == "image":
        mark = "x"
        col = ["y"]
        order = 0
    elif type == "source":
        mark = "s"
        col = ["k"]
        order = 0
    elif type == "mic":
        mark = "d"
        col = ["g"]
        order = 0
    elif type == None:
        mark = "."
        col = "k"
        order = 0

    ax.scatter(point.x, point.y, point.z, marker=mark, color=col[order])


def plot_line(ax, pointA, pointB):

    x, y, z = [pointA.x, pointB.x], [pointA.y, pointB.y], [pointA.z, pointB.z]
    ax.plot(x, y, z, color="black")
