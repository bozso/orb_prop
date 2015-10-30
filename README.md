# Orbital propagation using Octave/Matlab

## Goals

Comparing two numerical methods of satellite trajectory calculation: SGP4 and numerical integration of the equations of motion of a "low Earth orbit" satellite using the ode78 function. For that we use the implementation of the SGP4 model that can be downloaded here: (http://celestrak.com/software/tskelso-sw.asp). The numerical integrator is under development and can be accessed from this repository.  

## Input of the integrator

Prameter file should look like this:  
step day model  
orbital elements

## Currently working on:
* gravitational acceleration determined from numerical derivation of the potential  
* run_rk78.m reading inital conditions from orbit element sets

## History
* 2015-10-30: Implenetation of get_ics.m in run_rk78.m has begun, now working on adding 
functionality to the script that will read inital conditions from orbital element sets
* 2015-10-29: Added get_ics.m script that takes an orbital element set and extracts
the inital conditions from it (positin and velicity vector)
