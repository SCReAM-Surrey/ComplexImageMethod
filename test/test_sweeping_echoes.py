#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 27 10:39:38 2022

@author: od0014
"""
import numpy as np
import os
import pyroomacoustics as pra
from pyroomacoustics.transform import stft
from pyroomacoustics import metrics as met
import matplotlib.pyplot as plt
from scipy.io.wavfile import write
from scipy.fft import ifft, fft

from cim import SimpleGeometry as geom
from cim import ComplexImageSimulation as sim
from cim.utils import visualize_room as vis

plt.rcParams.update({"font.size": 8})

# whether to plot room and mic config.
plot_room = True
# whether to save plot
save = True
audio_path = "audio/"
fig_path = "figures/"

if save:
    if not os.path.exists(audio_path):
        os.mkdir(audio_path)
    if not os.path.exists(fig_path):
        os.mkdir(fig_path)

# =============================================================================
# helper function
# =============================================================================
def plot_spectrogram(
    rir,
    fs,
    fft_size,
    fft_hop,
    fft_zp,
    analysis_window,
    method,
    save=False,
    absorp=0,
):

    S = stft.analysis(rir, fft_size, fft_hop, win=analysis_window, zp_back=fft_zp)
    fig, (ax1, ax2) = plt.subplots(2, 1)

    if save:
        fig.set_size_inches(4.0, 3.0)

    ax1.imshow(
        pra.dB(S.T),
        extent=[0, len(rir), 0, fs / 2],
        vmin=-100,
        vmax=0,
        origin="lower",
        cmap="jet",
    )
    ax1.set_title("RIR generated with " + method)
    ax1.set_ylabel("Frequency")
    ax1.set_aspect("auto")

    # plot RIR
    ax2.plot(rir)
    ax2.set_xlabel("Num samples")
    ax2.set_ylabel("Amplitude")

    if save:
        plt.savefig(
            fig_path + "spec_" + method + "_absorp=" + str(absorp) + ".png", dpi=1000
        )

    plt.show()


def absorp_to_admittance(absorp):
    # find reflection coefficient
    reflect = np.sqrt(1 - absorp)
    # r = (Z/Za - 1/(Z/Za + 1)
    # Z = impedance, beta = admittance
    beta = (1.0 - reflect) / (1.0 + reflect)
    return beta


def admittance_to_absorp(beta):

    beta = np.real(beta)
    r = (1 - beta) / (1 + beta)
    absorp = 1 - np.power(r, 2)
    return absorp


# =============================================================================
# Run the example
# =============================================================================
# create an example with sweeping echo - from Enzo's paper
room_size = [4, 4, 4]
source_loc = [1, 2, 2]
center_x = room_size[0] / 2
center_y = room_size[1] / 2
center_z = room_size[2] / 2

# number of mics
M = 15
mic_loc = np.zeros((3, M), dtype=float)
gridSpacing = 0.25


# absorption coefficient
absorp = 0.5

# wall admittance
beta = absorp_to_admittance(absorp)

# sampling frequency
fs = 4000

# order for IM
ref_order = 10

# =============================================================================
# create uniform grid of mics
# =============================================================================
count = 0
Nmic_each_axis = int((M / 3 - 1) / 2)
for n in range(-Nmic_each_axis, Nmic_each_axis + 1):
    mic_loc[:, count] = [center_x + n * gridSpacing, center_y, center_z]
    mic_loc[:, count + 1] = [center_x, center_y + n * gridSpacing, center_z]
    mic_loc[:, count + 2] = [center_x, center_y, center_z + n * gridSpacing]
    count += 3

# which receiver position to investigate
which_receiver = 0

# =============================================================================
# create ComplexImageMethod objects
# =============================================================================

srcPos = geom.Point(source_loc[0], source_loc[1], source_loc[2])

room = geom.Room()
room.shape = geom.Cuboid(room_size[0], room_size[1], room_size[2])

receiverPos = list()

for i in range(M):
    receiverPos.append(geom.Point(mic_loc[0, i], mic_loc[1, i], mic_loc[2, i]))


# =============================================================================
# # draw the setup
# =============================================================================
# plot to visualize room - optional

if plot_room:
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

    for k in range(M):
        vis.plot_point(ax, receiverPos[k], "mic")

    ax.view_init(45, 110)
    if save:
        plt.savefig(fig_path +"test_sweep_setup.png", dpi=1000)
    plt.show()


# =============================================================================
# # normal ISM with sweeping echoes
# =============================================================================

# an absorption coefficient of 0 means purely reflective walls
room_ism = pra.ShoeBox(
    room_size, fs, materials=pra.Material(absorp), max_order=ref_order
)

room_ism.add_source(source_loc)

room_ism.add_microphone_array(pra.MicrophoneArray(mic_loc, fs))

room_ism.unset_air_absorption()

room_ism.compute_rir()

# plot spectrograms to check for sweeping echoes

fft_size = 64  # fft size for analysis
fft_hop = 16  # hop between analysis frame
fft_zp = 16  # zero padding
analysis_window = pra.hann(fft_size)


rir_ism = room_ism.rir[which_receiver][0]
Nfft = len(rir_ism)
ism_fft = fft(rir_ism)
freqs = np.linspace(0, fs / 2.0, int(Nfft / 2))

# ssf = met.sweeping_echo_measure(rir_ism,fs, t_min=0, t_max=0.5, fb=100)
# print('ISM sweeping echo measure' , ssf)

plot_spectrogram(
    rir_ism,
    fs,
    fft_size,
    fft_hop,
    fft_zp,
    analysis_window,
    "ISM",
    save,
    np.round(absorp, 3),
)

if save:
    write(
        audio_path + "ISM_absorp="
        + str(np.round(absorp, 3))
        + "_order="
        + str(ref_order)
        + ".wav",
        fs,
        rir_ism,
    )

# =============================================================================
# try the complex image method
# =============================================================================

# must be at least twice the number of samples in RIR

Nfreqs = 4000
#  speed of sound in air
c = 343
# wave numbers to evaluate on, equivalent to frequencies 0 to fs/2
wave_nums = np.linspace(0.01, np.pi * fs / c, Nfreqs)


# add wall impedance
nWalls = room.shape.nWalls
wallNames = ["floor", "ceiling", "left", "right", "front", "back"]


for k in range(nWalls):
    # constant admittance
    room.wallImpedance[wallNames[k]] = [beta * (1 + 1j) for i in range(Nfreqs)]


csim = sim.ComplexImageSimulation(room, srcPos, receiverPos, wave_nums, order=ref_order)
ref_pressure, total_pressure = csim.run()

pressure_at_receiver = np.array(total_pressure[which_receiver, :], dtype="complex")

rir = ifft(np.r_[pressure_at_receiver, np.flip(np.conj(pressure_at_receiver[1:]))])

# try to make the response minimum phase
rir_cism = np.real(rir[: len(rir_ism)])


# =============================================================================
# compare and make plots
# =============================================================================

# plot frequency spectrum
fig, ax = plt.subplots()
fig.set_size_inches(3.3, 2.5)
(ax1,) = ax.semilogx(
    wave_nums * (c / (2 * np.pi)), 20 * np.log10(np.abs(pressure_at_receiver))
)
(ax2,) = ax.semilogx(freqs, 20 * np.log10(np.abs(ism_fft[: int(Nfft / 2)])))
plt.xlim([10, fs / 2.0])
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude (dB)")
ax.legend([ax1, ax2], ["CISM", "ISM"])
if save:
    plt.savefig(fig_path + "freq_comp_absorp=" + str(absorp) + ".png", dpi=1000)
plt.show()

#  plot spectrograms
plot_spectrogram(
    rir_cism,
    fs,
    fft_size,
    fft_hop,
    fft_zp,
    analysis_window,
    "CISM",
    save,
    np.round(absorp, 3),
)

if save:
    write(
        audio_path + "CISM_absorp="
        + str(np.round(absorp, 3))
        + "_order="
        + str(ref_order)
        + ".wav",
        fs,
        rir_cism,
    )

# ssf_c = met.sweeping_echo_measure(rir_cism,fs, t_min=0, t_max=0.5, fb=100)
# print('CISM sweeping echo measure' , ssf_c)

# =============================================================================
# # the position of the direct path, as given by CSIM is correct. PRA is adding
# an additional delay
# =============================================================================
