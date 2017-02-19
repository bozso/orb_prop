% Initalizes necessary paths and contants

global tumin mu mu_si radiusearthkm a_earth xke C whichconst day2sec ...
sec2day eod deltaT_jday omega

omega = [0 0 1.15740740740741e-05];

% Acces to function scripts
addpath ('functions');
addpath ('vallado');
addpath('/usr/share/octave/packages/odepkg-0.8.5');


% Turn off annoying warning message
warning('off', 'all');

eod = load('eod_IAU2000_MJD.dat');
load deltaT_jday;

% WGS84
whichconst = 84;

% Necessary constants
C = [-0.484165143790815e-3, 0.95716107093473e-6, 0.53996586663899e-6];

% Getting WGS84 constants
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(whichconst);

% Changing into SI units
a_earth = radiusearthkm * 1000;
mu_si = mu * 1e9;

day2sec = 86400;
sec2day = 1 / day2sec;
