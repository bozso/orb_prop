% Initalizes necessary paths and contants

% Acces to function scripts
addpath ('/home/istvan/orb_prop/functions');
addpath ('/home/istvan/orb_prop/sgp4');

% Turn off annoying warning message
warning('off', 'Octave:possible-matlab-short-circuit-operator');

global tumin mu mu_si radiusearthkm a_earth xke C whichconst

% WGS84
whichconst = 84;

% Necessary constants
C = [-0.484165143790815e-3, 0.95716107093473e-6, 0.53996586663899e-6];

% Getting WGS84 constants
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(whichconst);

% Changing into SI units
a_earth = radiusearthkm * 1000;
mu_si = mu * 1e9;
