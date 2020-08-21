#!/usr/bin/env python3
# ref. https://numpy.org/doc/1.18/numpy-user.pdf

# Single Beam Sonar
#  Sonar Point-Scatter Model
# Contributors: Andreina Rascon, Derek Olson, Woeng-Sug Choi

#import cv2
from random import random
from math import sqrt, sin, cos, pi, log
import numpy as np
import matplotlib.pyplot as plt

# sinc function
def custom_sinc(x):
    return np.abs(np.sin(x)/x)

# Input Parameters
# The BlueView P900-45 is used as an example sonar for the purposes of
# demonstrating the point-scatter model

## Sonar Parameters
sonarFreq = 900E3 # Sonar frquency
bandwidth = 29.5e4 # [Hz]
freqResolution = 100e2
beamWidth = 0.1 # radians
nRays = 3

# Beam Parameters
nBeams = 7
np1_th = np.linspace(-beamWidth,beamWidth,nBeams)

## Input Parameters
soundSpeed = 1500.0 # [m/s]

# Target information (which should be obtained from ray gazebo plugin)
# 3 Targets, One dimension, no incident angle
np1_ray_distance = [1, 5, 1]
np1_ray_alpha = [0, 0, 0]
np1_ray_azimuthAngleWidth = [0.1, 0.1, 0.1]
np1_ray_elevationAngleWidth = [0.1, 0.1, 0.1]

# Surface Reflectivity
mu = 10e-4

# Frequency spectrum
fmin = sonarFreq - bandwidth/2*4
fmax = sonarFreq + bandwidth/2*4

# Sampling periods
max_T = max(np1_ray_distance)*2/soundSpeed
delta_f = 1/max_T
delta_t = 1/(fmax - fmin)

nFreq = round((fmax - fmin) / delta_f)
np1_freq = np.linspace(fmin,fmax,nFreq)
np1_time = np.linspace(0,max_T,nFreq)
NT = nFreq # dimensionless number of np1_time samples

# Transmit spectrum, frequency domain
np1_S_f = 1e11 * np.exp(-(np1_freq - sonarFreq)**2 * pi**2 / bandwidth**2)

# define wave vector
kw = 2 * pi * np1_freq / soundSpeed

# add attenuation
absorption = 0.0354 # [dB/m]
attenuation = absorption * log(10)/20
K = kw + 1j * attenuation

## Point Scattering model
# Calculate echo level using the point scatter model

np3_P_f = np.zeros((nBeams,nRays,nFreq), dtype=np.complex_)  # pre-allocate P(f)
np2_P_f_ray = np.zeros((nRays,nFreq), dtype=np.complex_) # pre-allocate ray P(f)
np2_P_t_modified = np.zeros((nBeams,nFreq))

# Beamforming Loop
for k in range(nBeams):

    # Ray Loop
    for i in range(nRays):

        xi_z = random() # [0,1) generate a random number, (Gaussian noise)
        xi_y = random() # [0,1) generate another random number, (Gaussian noise)

        alpha = np1_ray_alpha[i]
        distance = np1_ray_distance[i]
        azimuthAngleWidth = np1_ray_azimuthAngleWidth[i]
        elevationAngleWidth = np1_ray_elevationAngleWidth[i]
        amplitude = ((xi_z + 1j * xi_y) / sqrt(2)) * sqrt(
                         mu * cos(alpha)**2 * distance**2 \
                         * azimuthAngleWidth * elevationAngleWidth)

        # travel time
        travelTime = distance * 2 / soundSpeed


        # Summation of Echo returned from a signal
        np2_P_f_ray[i,:] = (np2_P_f_ray[i,:] \
                    + np1_S_f*amplitude*np.exp(1j*K*distance*2)/(distance**2))


    np3_P_f[k,i,:] = np2_P_f_ray[i,:]
    np2_P_beam = np.sum(np3_P_f, 1)

    np1_beamPatternCurrent = np.nan_to_num(custom_sinc(np1_th / beamWidth)**2,
                                           nan=1.0)

    # frequency domain to time domain
    print("np2_P_beam")
    print(np2_P_beam[:,10])
    np2_P_t = np.fft.ifft(np2_P_beam,axis=1)
    print("np2_P_t")
    print(np2_P_t[:,10])

    # p_t is nFreqxnBeams,
    # and np1_beamPatternCurrent is nBeamsx1
    np2_P_t_modified = np1_beamPatternCurrent.T[:,None] * np2_P_t

## Plotting
# Plot inverse fast fourier transform
plt.figure(figsize=(10,8))
plt.suptitle("%s rays, %d beams"%(nRays, nBeams))

# figure (1)
plt.subplot(2,2,1)
for i in range(nBeams):
    plt.plot(np1_time,np2_P_t[i,:], linewidth=0.5)
plt.xlabel('Time, [s]')
plt.ylabel('Pressure, [Pa]')

# # Plot Sound Pressure Level of Echo Level
np2_SPL = 20*np.log(np.abs(np2_P_t)) # sound pressure level, [dB]

# figure (2)
plt.subplot(2,2,2)
for i in range(nBeams):
    plt.plot(np1_time,np2_SPL[i,:], linewidth=0.5)
plt.xlabel('Time, [s]')
plt.ylabel('Sound Pressure Level, [dB]')

# # Plotting Beam Formed Data
# sound pressure level, [dB]
np2_SPL_modified = 20*np.log(np.abs(np2_P_t_modified))
print("np2_SPL_modified")
print(np2_SPL_modified[:,10])

# figure (3)
plt.subplot(2,2,3)
for i in range(nBeams):
    plt.plot(np1_time,np2_SPL_modified[i,:], linewidth=0.5)
plt.xlabel('Time, [s]')
plt.ylabel('Sound Pressure Level modified, [dB]')

# figure (4)
plt.subplot(2,2,4)
plt.imshow(np2_SPL_modified.T, aspect="auto")
plt.show()

