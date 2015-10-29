# Orbital propagation using Octave/Matlab

## Goals

Comparing two numerical methods of satellite trajectory calculation: SGP4 and numerical integration of the equations of motion of a "low Earth orbit" satellite using the ode78 function. For that we use the implementation of the SGP4 model that can be downloaded here: (http://celestrak.com/software/tskelso-sw.asp). The numerical integrator is under development and can be accessed from this repository.  

## Input of the integrator

Prameter file should look like this:  
t0      t0 time  
step    # of seconds in a section  
day     # of days to propagate  
name    name of output files  
inital conditions x, y, z, vx, vy, vz separated by spaces

## Currently working on:
* gravitational acceleration determined from numerical derivation of the potential  

## History
* 2015-10-29: Added gt_ics.m script that takes an orbital element set and extracts
the inital conditions from it (positin and velicity vector)
