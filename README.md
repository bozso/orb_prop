# Orbital propagation using Octave/Matlab

## Goals

Comparing two numerical methods of satellite trajectory calculation: SGP4 and numerical integration of the equations of motion of a "low Earth orbit" satellite using the ode78 function. For that we use the implementation of the SGP4 model that can be downloaded here: (http://celestrak.com/software/tskelso-sw.asp). The numerical integrator is under development and can be accessed from this repository.  

## Input of the integrator

Prameter file should look like this:  
step day model  
orbital elements

## Currently working on:
* comparing sgp4 and Runge-Kutta propagation without friction
* integrating air friction into the Runge-Kutta model
