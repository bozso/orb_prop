
%~ Processing parameters and initial conditions
[step day model satrec] = process_infile('/home/istvan/orb_prop/synthetic.tle', whichconst);

switch (model)
	case 'point_mass'
		f_handle = @point_mass;
		nmodel = 1;
	case 'zonal'
		f_handle = @zonal;
		nmodel = 2;
	otherwise
		disp('Unrecognized model option!');
		return
end

vopt = odeset('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', 'InitialStep', 100, 'MaxStep', 100);

%~ Getting initial conditions
[satrec, ro ,vo] = sgp4(satrec,  0.0);

%~ Changing into SI units
rv_0 = [ro, vo] * 1e3;

%~ plotname is satellite number
plotname = ['output/', num2str(satrec.satnum)];


t0 = 0.0;

if (nmodel == 1)
	%~ Calculation  of orbital element set from initial condtions
	[a, eccentricity, inclination, omega, sinw, cosw, n, tau] = ...
	orbital_elements(t0, rv_0);
	
	%~ Calculating Rot matrix that transforms coordinates from the orbital plane to the 
	%~ reference plane.
	Rot = calc_rotmatrix(sinw, cosw, inclination, omega);
end



%~ Propagation with ode78
[t, r_ode] = ode78(f_handle, [t0, 86400 * day], rv_0, vopt);

timexyzv = [t, r_ode];

save('-ascii', [plotname, '.dat'], 'timexyzv');
clear timexyzv;
	
%~ Getting x, y, z coordinates
xyz = [r_ode(:,1), r_ode(:,2), r_ode(:,3)];

if (nmodel == 1)
	%~ Calculation of the theoretical coordinates in the orbital plane
	ellipse = calc_ellipse (t, tau, eccentricity, a, n);

	ellipse = ellipse * Rot;

	save('-ascii', [plotname, '_ellipse.dat'], 'ellipse');

	%~ Calculating differences between theoretical and numerical coordinates
	delta = sqrt(sum((ellipse - xyz)'.^2));
	
	delta = [t, delta'];
	save('-ascii', [plotname, '_delta.dat'], 'delta');

	%~ Plotting deltas
	h = figure(2);
	clf
	hold on
	plot(delta(:,1) / 86400, delta(:,2));
	print(h,'-dpng','-color', [plotname, '_delta' , '.png']);
	hold off
end % if point_mass

%~ Plotting orbit
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
