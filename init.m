%~ Acces to function scripts
addpath ('/home/istvan/orb_prop/functions');
addpath ('/home/istvan/orb_prop/sgp4');

warning('off', 'Octave:possible-matlab-short-circuit-operator');

%~ WGS84
whichconst = 84;

%~ Necessary constants
global tumin mu mu_si radiusearthkm a_earth xke C
C = [-0.484165143790815e-3, 0.95716107093473e-6, 0.53996586663899e-6];

%~ Getting WGS84 constants
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(whichconst);

%~ Changing into SI units
a_earth = radiusearthkm * 1000;
mu_si = mu * 1e9;
