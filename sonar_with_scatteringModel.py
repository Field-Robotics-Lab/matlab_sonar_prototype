#!/usr/bin/env python3
# ref. https://numpy.org/doc/1.18/numpy-user.pdf

# Single Beam Sonar
#  Sonar Point-Scatter Model
# Contributors: Andreina Rascon, Derek Olson, Woeng-Sug Choi

from random import random
from math import sqrt, sin, cos, pi, log
import numpy as np
import matplotlib.pyplot as plt
import time

time_start = time.perf_counter() 

# unnormalized sinc function
def unnormalized_sinc(t):
    return np.nan_to_num(np.sin(t)/t, nan=1.0)

## physics constants
soundSpeed = 1500.0 # [m/s]
mu = 10e-4          # Surface Reflectivity

# Input Parameters
# The BlueView P900-45 is used as an example sonar for the purposes of
# demonstrating the point-scatter model

# Sonar properties
sonarFreq = 900E3 # Sonar frquency
bandwidth = 29.5e4 # [Hz]
freqResolution = 100e2
fmin = sonarFreq - bandwidth/2*4 # Calculated requency spectrum
fmax = sonarFreq + bandwidth/2*4 # Calculated requency spectrum

# Sonar sensor properties
nBeams = 1
beam_elevationAngle = 0.0175 # Beam looking down in elevation direction
beam_azimuthAngle = 0.0 # Beam at center line in azimuth direction
beam_elevationAngleWidth = 0.1 # radians
beam_azimuthAngleWidth = 0.1 # radians
ray_nElevationRays = 4
ray_nAzimuthRays = 3

# calculated Sonar sensor properties
ray_elevationAnglesf1 = beam_elevationAngle + np.linspace(
                 -beam_elevationAngleWidth / 2, beam_elevationAngleWidth / 2,
                 ray_nElevationRays)
ray_azimuthAnglesf1 = beam_azimuthAngle + np.linspace(
                 -beam_azimuthAngleWidth / 2, beam_azimuthAngleWidth / 2,
                 ray_nAzimuthRays)
ray_elevationAngleWidth = beam_elevationAngleWidth/(ray_nElevationRays - 1)
ray_azimuthAngleWidth = beam_azimuthAngleWidth/(ray_nAzimuthRays - 1)

# Sonar inputs
ray_distancef2 = np.array([[15,5,10], [2,10,10], [15,15,15], [4,2,3]])
ray_alphaf2 = np.array([[0,0,0], [0,0,0], [0,0,0], [0,0,0]]) # no incident angle

# calculated sampling periods
_max_T = np.amax(ray_distancef2)*2/soundSpeed
_delta_f = 1/_max_T
# _delta_t = 1/(fmax - fmin)
nFreq = round((fmax - fmin) / _delta_f)
_freq1f = np.linspace(fmin,fmax,nFreq)
time1f = np.linspace(0,_max_T,nFreq) # for diagnostics plot

# calculated physics
_absorption = 0.0354 # [dB/m]
_attenuation = _absorption*log(10)/20
_kw3f = 2*pi*_freq1f/soundSpeed # wave vector
K3f = _kw3f + 1j*_attenuation # attenuation constant K3f

# Transmit spectrum, frequency domain
S_f1f = 1e11 * np.exp(-(_freq1f - sonarFreq)**2 * pi**2 / bandwidth**2)

# Point Scattering model
# Echo level using the point scatter model for P(f) and P(t) for beams
P_ray_f3c = np.zeros((ray_nElevationRays,ray_nAzimuthRays,nFreq),
                      dtype=np.complex_)
# P_ray_t3f = np.zeros((ray_nElevationRays,ray_nAzimuthRays,nFreq),
#                      dtype=np.complex_)
azimuthBeamPattern2f = np.zeros((ray_nElevationRays,ray_nAzimuthRays))
elevationBeamPattern2f = np.zeros((ray_nElevationRays,ray_nAzimuthRays))
for k in range(ray_nElevationRays):
    for i in range(ray_nAzimuthRays):
        azimuthBeamPattern2f[k,i] = (np.abs(unnormalized_sinc(pi * 0.884
               / ray_azimuthAngleWidth * sin(ray_azimuthAnglesf1[i]))))**2
        elevationBeamPattern2f[k,i] = (np.abs(unnormalized_sinc(pi * 0.884
               / ray_elevationAngleWidth * sin(ray_elevationAnglesf1[k]))))**2

for k in range(ray_nElevationRays):
    for i in range(ray_nAzimuthRays):
        xi_z = random()   # generate a random number, (Gaussian noise)
        xi_y = random()   # generate another random number, (Gaussian noise)

        # angle between ray vector and object normal vector, [rad]
        alpha = ray_alphaf2[k,i]

        distance = ray_distancef2[k,i]
        amplitude = (((xi_z + 1j * xi_y)
                     / sqrt(2))
                     * (sqrt(mu * cos(alpha)**2 * distance**2
                             * ray_azimuthAngleWidth
                             * ray_elevationAngleWidth))
                     * azimuthBeamPattern2f[k,i]
                     * elevationBeamPattern2f[k,i])

        # Summation of Echo returned from a signal (frequency domain)
        for m in range(nFreq):
            P_ray_f3c[k,i,m] = P_ray_f3c[k,i,m] + S_f1f[m] * amplitude \
                       * np.exp(-1j * K3f[m] * distance * 2) / (distance**2)

# Summation for a beam
P_beamf1 = np.sum(np.sum(P_ray_f3c, 1), 0)
P_beam_tf1 = np.fft.ifft(P_beamf1)

# Calculate computation time
time_elapsed = (time.perf_counter()  - time_start)
print(time_elapsed)

## Plots
plt.figure(figsize=(10,8))
plt.suptitle("%s elevation rays, %d azimuth rays"%(ray_nElevationRays,
                                                   ray_nAzimuthRays))

# inverse fast fourier transform
# figure (1)
plt.subplot(2,1,1)
plt.grid(True)
plt.plot(time1f, P_beam_tf1, linewidth=0.5)
plt.xlabel('Time, [s]')
plt.ylabel('Pressure, [Pa]')

# Sound Pressure Level of Echo Level
# figure (2)
SPLf1 = 20 * np.log(np.abs(P_beam_tf1)) # sound pressure level, [dB]
plt.subplot(2,1,2)
plt.grid(True)
plt.plot(time1f, SPLf1, linewidth=0.5)
plt.xlabel('Time, [s]')
plt.ylabel('Sound Pressure Level, [Pa]')

plt.show()

