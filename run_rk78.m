% Acces to function scripts
addpath ('/home/istvan/orb_prop/rk_orbit/functions');
addpath ('/home/istvan/orb_prop/sgp4');

% Necessary constants
global tumin mu a xke j
[tumin, mu, a, xke, j2, j3, j4, j3oj2] = getgravc(84);
j = [j2, j3, j4];
a = a * 1000;
mu = mu * 1e9;


% File that contains propagation parameters and initial conditions
fin = fopen ('ics.dat');

% Processing parameters and initial conditions
for iii = 1:5
	str = fgetl(fin);
	switch ( str(1:findstr(str, ' ') - 1)  )
		case 't0'
			t0 = str2num (str(rindex(str, ' '):end));
		case 'step'
			step = str2num (str(rindex(str, ' '):end));
		case 'day'
			day = str2num (str(rindex(str, ' '):end));
		case 'name'
			name = str(rindex(str, ' ') + 1:end);
		case 'model'
			model = str(rindex(str, ' ') + 1:end);
		otherwise
			disp('Warning: Unrecognized option.')
	end
end

str = fgetl(fin);
rv_0 = strread(str);

fclose(fin);

plotname = ['output/', name];

vopt = odeset ('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', 'InitialStep', 125, 'MaxStep', 200);

% Calculation  of orbital element set from initial condtions
[a, eccentricity, inclination, omega, sinw, cosw, n, tau] = OrbitalElementSet (rv_0, t0, GRAVITY_CONST * EARTH_MASS);

% Calculating Rot matrix that transforms coordinates from the orbital plane to the 
% reference plane.
Rot = CalcRotMtx (sinw, cosw, inclination, omega);

% How many 'step' seconds long sections we need
cicles = ceil((86400 * day ) / step);

% Propagation with ode78
timexyzv = rk78 (rv_0, step, cicles, model, vopt);

% Getting x, y, z coordinates
xyz = [timexyzv(:,2), timexyzv(:,3), timexyzv(:,4)];

% Getting time points
time = timexyzv(:, 1);

save('-ascii', [plotname, '.dat'], 'timexyzv');

clear timexyzv;

% Calculation of the theoretical coordinates in the orbital plane
ellipse = CalcEllipse (time, tau, eccentricity, a, n);

ellipse = ellipse * Rot;

save ('-ascii', [plotname, '_ellipse.dat'], 'ellipse');

% Calculating differences between theoretical and numerical coordinates
delta = ellipse - xyz;

delta = sqrt(sum(delta'.^2));

delta = [time, delta'];

% Plotting orbit
h=figure(1);
hold on;
plot3 (xyz(:, 1) / 1e6, xyz(:, 2) / 1e6, xyz(:, 3) / 1e6, 'ob;RK4;', 'linewidth', 1, ...
ellipse(:,1) / 1e6, ellipse(:,2) / 1e6, ellipse(:,3) / 1e6, '*r;Ellipse;', 'linewidth', 1, ...
0, 0, 0, 'ok;Earth;');
view(3);
xlabel('x [1000 km]');
ylabel('y [1000 km]');
zlabel('z [1000 km]');
legend ('location', 'northeastoutside');
print(h,'-dpng','-color', [plotname, '.png']);
hold off

save('-ascii', [plotname, '_delta.dat'], 'delta');

% Plotting deltas
h = figure(2);
clf
hold on
plot (delta(:,1) / 86400, delta(:,2));
print(h,'-dpng','-color', [plotname, '_delta' , '.png']);
hold off

