# Orbital propagation using Octave/Matlab

## Goals

Comparing two numerical methods of satellite trajectory calculation:
SGP4 and numerical integration of the equations of motion of a "low Earth orbit"
satellite using the ode78 function.
For that we use the implementation of the SGP4 model that can be downloaded here:
(http://celestrak.com/software/tskelso-sw.asp).
The numerical integrator is under development and can be accessed from this repository.
The codes use the SGP library
for setting the initial parameters of integration.

## Input of the integrator

Parameter file should look like this:  
step day model  
orbital elements  

For example:  
100 1 zonal  
1 11416U 79057A   80003.44366214  .00000727  00000-0  33454-3 0   878
2 11416  98.7309  35.7226 0013335  92.0280 268.2428 14.22474848 27074
  
Stepsize = 100 seconds, integration time = 1 day. More on the two-line element set:
(https://celestrak.com/columns/v04n03/)  

## Currently working on:
* comparing sgp4 and Runge-Kutta propagation without friction
* integrating air friction into the Runge-Kutta model
