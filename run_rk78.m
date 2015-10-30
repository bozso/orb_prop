%~ Acces to function scripts
addpath ('/home/istvan/orb_prop/rk_orbit/functions');
addpath ('/home/istvan/orb_prop/sgp4');

%~ Necessary constants
global tumin mu mu_si radiusearthkm a xke j j2 j3 j4 j3oj2

%~ Getting WGS84 constants
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(84);
j = [j2, j3, j4];

%~ Changing into SI units
a = radiusearthkm * 1000;
mu_si = mu * 1e9;

%~ File that contains propagation parameters and initial conditions
infile = fopen (argv(){1});

%~ Processing parameters and initial conditions

param = strread(fgetl(infile));

step = param(1);
day = param(2);

%~ Processing oribtal element sets
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
	printf('Reading satellite data: %d\r', iii);
	fflush(stdout);
end

fclose(infile);



vopt = odeset ('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', 'InitialStep', 125, 'MaxStep', 200);

for jjj = 1:(iii - 1)
	
	[satrec(jjj), ro ,vo] = sgp4 (satrec(jjj),  0.0);
	
	rv_0 = [ro * 1e3, vo * 1e3];
	
	plotname = ['output/', num2str(satrec(jjj).satnum)];
	
	if (model == 'point_mass')
		%~ Calculation  of orbital element set from initial condtions
		[a, eccentricity, inclination, omega, sinw, cosw, n, tau] = OrbitalElementSet (rv_0, t0, GRAVITY_CONST * EARTH_MASS);
	end
	
	%~ Calculating Rot matrix that transforms coordinates from the orbital plane to the 
	%~ reference plane.
	Rot = CalcRotMtx (sinw, cosw, inclination, omega);

	%~ How many 'step' seconds long sections we need
	cicles = ceil((86400 * day ) / step);

	%~ Propagation with ode78
	timexyzv = rk78 (rv_0, step, cicles, model, vopt);

	%~ Getting x, y, z coordinates
	xyz = [timexyzv(:,2), timexyzv(:,3), timexyzv(:,4)];

	%~ Getting time points
	time = timexyzv(:, 1);

	save('-ascii', [plotname, '.dat'], 'timexyzv');

	clear timexyzv;

	if (model == 'point_mass')
		%~ Calculation of the theoretical coordinates in the orbital plane
		ellipse = CalcEllipse (time, tau, eccentricity, a, n);
	
		ellipse = ellipse * Rot;

		save ('-ascii', [plotname, '_ellipse.dat'], 'ellipse');

		%~ Calculating differences between theoretical and numerical coordinates
		delta = ellipse - xyz;

		delta = sqrt(sum(delta'.^2));

		delta = [time, delta'];
		save('-ascii', [plotname, '_delta.dat'], 'delta');

		%~ Plotting deltas
		h = figure(2);
		clf
		hold on
		plot (delta(:,1) / 86400, delta(:,2));
		print(h,'-dpng','-color', [plotname, '_delta' , '.png']);
		hold off
	end
	
	%~ Plotting orbit
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
end %~ for - satrec
