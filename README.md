# Orbital propagation using Octave/Matlab

## Goals

Comparing two numerical methods of satellite trajectory calculation:
SGP4 and numerical integration of the equations of motion ofa "low Earth orbit" satellite using the ode78 function.
For that we use the implementation of the SGP4 model that can be downloaded here:
(http://celestrak.com/software/tskelso-sw.asp).
The numerical integrator is under development and can be accessed from this repository. The codes use the SGP librabry
for setting the initial parameters of integration.

## Input of the integrator

Prameter file should look like this:  
step day model  
orbital elements  

For example:  
100 1 zonal  
1 04793U 70106A   16038.50388155 +.00000000 +00000-0 +00000-0 0 00001  
2 04793 000.0000 000.0000 0032224 000.0000 000.0000 12.53975678000007  
  
Stepsize = 100 seconds, integration time = 1 day. More on the two-line element set:
(https://celestrak.com/columns/v04n03/)  

## Currently working on:
* comparing sgp4 and Runge-Kutta propagation without friction
* integrating air friction into the Runge-Kutta model
