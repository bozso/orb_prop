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
* a document that explains how the gravitational acceleration can be derived from the series expansion of the gravitational potential
* implementation of the gravitational acceleration, that is derived from the series expansion, in the script files 
