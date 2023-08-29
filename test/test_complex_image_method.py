#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 11:42:35 2022

@author: od0014
"""

import numpy as np
import os
from cim import ComplexImageSimulation as sim
from cim import SimpleGeometry as geom
from cim.utils import visualize_room as vis

from matplotlib import cm
import matplotlib.pyplot as plt


plt.rcParams.update({"font.size": 8})
# whether to plot the room and mic setup
plot_room = True
# whether to save the plots
save = True
# constant or frequency dependent admittance
admit = "soft"
# complex IM or normal IM
img_type = "complex"

# cuboid room dimensions
L = 4.0
W = 5.0
H = 4.5

# number of wave numbers
Nx = 75
# number of microphone positions
Ny = 200

srcPos = geom.Point(0.5, 2.3, 1.2)
sourceAngle = np.arctan(srcPos.y / srcPos.x)
wave_nums = np.linspace(1, 30, Nx)

# receivers are in an arc around xz plane
r = 2
theta = np.linspace(0, np.pi / 2, Ny)
rx = r * np.cos(theta)
ry = L / 2 * np.ones(Ny)
rz = r * np.sin(theta)
receiverPos = list()
srcRecAngle = np.zeros(Ny)

for i in range(Ny):
    curReceiverPos = geom.Point(rx[i], ry[i], rz[i])
    srcRecAngle[i] = srcPos.getAngle(curReceiverPos)
    receiverPos.append(curReceiverPos)

######################################
# run the simulation
orders = [1, 2, 5]

for order in orders:

    room = geom.Room()
    room.shape = geom.Cuboid(L, W, H)

    # add wall impedance
    # loop over 6 walls
    nWalls = room.shape.nWalls
    wallNames = ["floor", "ceiling", "left", "right", "front", "back"]

    for k in range(nWalls):
        if admit == "const":
            # constant admittance
            room.wallImpedance[wallNames[k]] = [0.5 * (1 + 1j) for i in range(Nx)]
        elif admit == "rigid":
            room.wallImpedance[wallNames[k]] = [0 for i in range(Nx)]
        else:
            # frequency-dependent admittance
            room.wallImpedance[wallNames[k]] = [
                (0.5 + 0.5 * 1j) * (1.0 / (1 + np.exp(-0.1 * wave_nums[i])))
                for i in range(Nx)
            ]

    csim = sim.ComplexImageSimulation(room, srcPos, receiverPos, wave_nums, order, img_type=img_type)
    ref_pressure, total_pressure = csim.run()

    #####################################
    # draw the setup
    if plot_room and order == 1:
        # plot to visualize room - optional
        fig = plt.figure()
        fig.set_size_inches(3.3, 2.5)
        ax = fig.add_subplot(projection="3d")
        vis.plot_room(
            ax,
            [room.shape.x / 2, room.shape.y / 2, room.shape.z / 2],
            room.shape.x,
            room.shape.y,
            room.shape.z,
        )
        vis.plot_point(ax, srcPos, "source")

        for k in range(Ny):
            vis.plot_point(ax, receiverPos[k], "mic")

        ax.text(2, 3, 0, r"$\theta = 0^{\circ}$")
        ax.text(1, 1, 1.5, r"$\theta = 90^{\circ}$")
        ax.view_init(45, 110)
        if save:
            if not os.path.exists('figures/'):
                os.mkdir('figures/')
            plt.savefig("figures/test1_setup.png", dpi=1000)
        plt.show()

    # surface plot

    [K, Theta] = np.meshgrid(wave_nums, theta)
    fig = plt.figure()
    fig.set_size_inches(3.3, 2.5)
    ax = fig.add_subplot(projection="3d")
    surf = ax.plot_surface(
        K,
        np.rad2deg(Theta),
        20 * np.log10(np.abs(total_pressure)),
        cmap=cm.jet,
        linewidth=0.1,
    )
    # if order == 2:
    #     fig.colorbar(surf, shrink=0.5, aspect=5, label="Pressure (dB)",
    #         orientation="vertical", labelpad=0.01)
    ax.set_xlabel("Wave number")
    ax.set_ylabel("Angle of receiver from origin")
    ax.set_zticks([])
    ax.view_init(270, -90)
    if save:
        plt.savefig(f"figures/admit={admit}_order={order}_surf.eps", format="eps")
    plt.show()

    # polar plot
    idx = [4, 10, 30]
    fig = plt.figure()
    for j in range(len(idx)):
        polar_p = total_pressure[:, idx[j]]
        ax = fig.add_subplot(len(idx), 1, j + 1, polar=True)
        ax.plot(theta, np.abs(polar_p))
        ax.set_thetamin(0)
        ax.set_thetamax(90)
        ax.set_title(
            "k = " + str(round(wave_nums[idx[j]], 3)),
            fontdict={"fontsize": 8, "fontweight": "light"},
        )
    # set the spacing between subplots
    plt.subplots_adjust(
        left=0.1, bottom=0.1, right=0.9, top=0.8, wspace=0.4, hspace=1.0
    )
    plt.suptitle(f"{admit} wall, order={order}")

    if save:
        plt.savefig(f"figures/admit={admit}_order={order}_polar.eps", format="eps")
    plt.show()
