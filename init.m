% Initalizes necessary paths and contants

% Acces to function scripts
addpath ('/home/istvan/orb_prop/functions');
addpath ('/home/istvan/orb_prop/sgp4');

% Turn off annoying warning message
warning('off', 'Octave:possible-matlab-short-circuit-operator');

global tumin mu mu_si radiusearthkm a_earth xke C whichconst day2sec ...
sec2day eod deltaT_jday jday_low jday_high

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

% loading eod
load('/home/istvan/orb_prop/eod_jday.mat');
% loading deltaT_jday
load('/home/istvan/orb_prop/deltaT_jday.mat');

jday_low = 2448622.50000000;
jday_high= 2457570.50000000;
