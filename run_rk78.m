%Acces to function scripts
addpath ('/home/istvan/orb_prop/rk_orbit/functions');
addpath ('/home/istvan/orb_prop/sgp4');
warning('off', 'Octave:possible-matlab-short-circuit-operator');

whichconst = 84;

% Necessary constants
global tumin mu mu_si radiusearthkm a_earth xke C
C = [-0.484165143790815e-3, 0.95716107093473e-6, 0.53996586663899e-6];

% Getting WGS84 constants
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(whichconst);

% Changing into SI units
a_earth = radiusearthkm * 1000;
mu_si = mu * 1e9;
%mu_si = 3.9858874e14;

% File that contains propagation parameters and initial conditions
infile = fopen('/home/istvan/orb_prop/fengyun_1c.tle');

% Processing parameters and initial conditions

str = fgetl(infile);

space_index = findstr(str, ' ');

if (length(space_index) > 2)
	disp('Error: Too many options given in the first line!')
	return
end

step = str2num(str(1:space_index(1) - 1));
day = str2num(str( (space_index(1) + 1) : (space_index(2) - 1) ));
model = str ( (space_index(2) + 1):end );

% Processing oribtal element sets

iii = 1;

% disp('Reading satellite data...')
while (~feof(infile))
    longstr1 = fgets(infile, 130);
    while ( (longstr1(1) == '#') && (feof(infile) == 0) )
        longstr1 = fgets(infile, 130);
    end

    if (feof(infile) == 0)
        longstr2 = fgets(infile, 130);
    end
        
	satrec(iii) = process(whichconst, longstr1, longstr2);
	iii = iii + 1;
end

fclose(infile);

vopt = odeset('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', 'InitialStep', 100, 'MaxStep', 100);


for jjj = 1:(iii - 1)
	% Getting initial conditions
	[satrec(jjj), ro ,vo] = sgp4(satrec(jjj),  0.0);
	
	% Changing into SI units
	rv_0 = [ro, vo] * 1e3;
	
	% plotname is satellite number
	plotname = ['output/', num2str(satrec(jjj).satnum)];
	
	% t0 = satrec(jjj).jdsatepoch * 86400;
	t0 = 0.0;
	
	%if (model == 'point_mass')
		% Calculation  of orbital element set from initial condtions
		%[a, eccentricity, inclination, omega, sinw, cosw, n, tau] = ...
		%orbital_elements(t0, rv_0);
		
		% Calculating Rot matrix that transforms coordinates from the orbital plane to the 
		% reference plane.
		%Rot = calc_rotmatrix(sinw, cosw, inclination, omega);
	%end
	
	switch (model)
		case 'point_mass'
			f_handle = @point_mass;
		case 'zonal'
			f_handle = @zonal;
		otherwise
			printf('Unrecognized model option!\n');
			return
	end

	% Propagation with ode78
	[t, r_int] = ode78(f_handle, [t0, 86400 * day], rv_0, vopt);
	
	timexyzv = [t, r_int];
	save('-ascii', [plotname, '.dat'], 'timexyzv');
	clear timexyzv;
		
	% Getting x, y, z coordinates
	xyz = [r_int(:,1), r_int(:,2), r_int(:,3)];
	
	% Getting time points
	time = t;

	%if (model == 'point_mass')
		% Calculation of the theoretical coordinates in the orbital plane
		%ellipse = calc_ellipse (time, tau, eccentricity, a, n);
	
		%ellipse = ellipse * Rot;

		%save('-ascii', [plotname, '_ellipse.dat'], 'ellipse');

		% Calculating differences between theoretical and numerical coordinates
		%delta = sqrt(sum((ellipse - xyz)'.^2));
		
		%delta = [time, delta'];
		%save('-ascii', [plotname, '_delta.dat'], 'delta');

		% Plotting deltas
		%h = figure(2);
		%clf
		%hold on
		%plot(delta(:,1) / 86400, delta(:,2));
		%print(h,'-dpng','-color', [plotname, '_delta' , '.png']);
		%hold off
	%end % if point_mass
	
	% Plotting orbit
	h=figure(1);
	hold on;
	plot3(xyz(:, 1) / 1e6, xyz(:, 2) / 1e6, xyz(:, 3) / 1e6, 'ob;RK4;', 'linewidth', 1, ...
	0, 0, 0, 'ok;Earth;');
	view(3);
	xlabel('x [1000 km]');
	ylabel('y [1000 km]');
	zlabel('z [1000 km]');
	legend('location', 'northeastoutside');
	print(h,'-dpng','-color', [plotname, '.png']);
	hold off
end % for satrec
