#!/usr/bin/env python3
# Single Beam Sonar

from random import random
from math import sqrt, sin, cos, pi, log
import numpy as np
import matplotlib.pyplot as plt
import time

time_start = time.perf_counter() 
print(time.perf_counter())

# unnormalized sinc function
def unnormalized_sinc(t):
    return np.nan_to_num(np.sin(t)/t, nan=1.0)

## physics constants
soundSpeed = 1500.0 # [m/s]
mu = 10e-4          # Surface Reflectivity
_absorption = 0.0354 # [dB/m]
sourceLevel = 150 # Source level (not critical since no absolute ambient noise level considered) 
noiseAmp = 0.01

# Sonar inputs
ray_distancef2 = np.array([[15,5,10], [2,10,10], [15,15,15], [4,2,3]])
ray_alphaf2 = np.array([[0,0,0], [0,0,0], [0,0,0], [0,0,0]]) # no incident angle

# Input Parameters
# The BlueView P900-45 is used as an example sonar
# Sonar properties
sonarFreq = 900E3 # Sonar frquency
bandwidth = 29.5e4 # [Hz]
maxRange = np.amax(ray_distancef2) # Maximum unambiguous range
rangeRes = 0.025 # Range resolution (1 inch = 0.0254 meter)

# Transmitted Waveform properties
prf = soundSpeed/(2.0*maxRange);            # Pulse repetition frequency
pulse_width = 2.0*rangeRes/soundSpeed;      # Pulse width
pulse_bw = 1.0/pulse_width;                 # Pulse bandwidth
fs = 2.0*pulse_bw;                          # Sampling rate

# Generate rectangular pulse and time domain
Domain_legnth = int(fs/prf)
Pulse = np.zeros(Domain_legnth); 
for i in range(int(fs*pulse_width)):
    Pulse[i-1] = 1
time1f = np.linspace(0,Domain_legnth-1,Domain_legnth)/fs

# Random Noise 
for i in range(Domain_legnth):
    Pulse[i] = Pulse[i] + noiseAmp*random()

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


azimuthBeamPattern2f = np.zeros((ray_nElevationRays,ray_nAzimuthRays))
elevationBeamPattern2f = np.zeros((ray_nElevationRays,ray_nAzimuthRays))
for k in range(ray_nElevationRays):
    for i in range(ray_nAzimuthRays):
        azimuthBeamPattern2f[k,i] = (np.abs(unnormalized_sinc(pi * 0.884
               / ray_azimuthAngleWidth * sin(ray_azimuthAnglesf1[i]))))**2
        elevationBeamPattern2f[k,i] = (np.abs(unnormalized_sinc(pi * 0.884
               / ray_elevationAngleWidth * sin(ray_elevationAnglesf1[k]))))**2

Signal = np.ones((ray_nElevationRays,ray_nAzimuthRays,Domain_legnth))
for k in range(ray_nElevationRays):
    for i in range(ray_nAzimuthRays):
        # angle between ray vector and object normal vector, [rad]
        alpha = ray_alphaf2[k,i]
        # distance 
        distance = ray_distancef2[k,i]
        # beamPattern
        beamPattern = azimuthBeamPattern2f[k,i] * elevationBeamPattern2f[k,i]

        #-- Source level
        raySL = sourceLevel + 10*np.log10(beamPattern)

        #-- Transmission path loss -- TL = 20log10(r) + alpha*r
        rayTL = 20*np.log10(distance) + _absorption*distance

        #-- Target Strength -- simple cosine model
        rayTS = 10*np.log10(sqrt(mu * cos(alpha)**2 * distance**2 
                             * ray_azimuthAngleWidth * ray_elevationAngleWidth))
        
        #-- Summation for final power level--
        rayEL = raySL - 2*rayTL + rayTS

        #-- Signal Calculation --
        rayRetTime = distance*2/soundSpeed
        Signal[k,i] = np.roll(Pulse*rayEL,int(rayRetTime*fs))
        

# Summation for a beam
SPLf1 = np.sum(np.sum(Signal, 1), 0) 
SPLf1[0] = 0; SPLf1[len(SPLf1)-1] = 0
P_beam_tf1 = 10**(SPLf1/20)

# Calculate computation time
time_elapsed = (time.perf_counter()  - time_start)
print(time.perf_counter())
print(time_elapsed)

## Plots
plt.figure(figsize=(10,8))
plt.suptitle("%s elevation rays, %d azimuth rays"%(ray_nElevationRays,
                                                   ray_nAzimuthRays))

# figure (1)
plt.subplot(2,1,1)
plt.grid(True)
plt.plot(time1f, P_beam_tf1, linewidth=0.5)
plt.xlabel('Time, [s]')
plt.ylabel('abs(Pressure), [Pa]')

# Sound Pressure Level of Echo Level
# figure (2)
plt.subplot(2,1,2)
plt.grid(True)
plt.plot(time1f, SPLf1, linewidth=0.5)
plt.xlabel('Time, [s]')
plt.ylabel('Sound Pressure Level, [dB]')

plt.show()

