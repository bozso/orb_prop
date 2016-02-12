function [] = rk_propagate(infile_path)

	%~ Processing parameters and initial conditions
	[step day model satrec] = process_infile(infile_path, whichconst);

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

	%~ Propagation with ode78
	[t, r_ode] = ode78(f_handle, [t0, 86400 * day], rv_0, vopt);

	timexyzv = [t, r_ode];

	save('-ascii', [plotname, '.dat'], 'timexyzv');
	clear timexyzv;
		
	%~ Getting x, y, z coordinates
	xyz = [r_ode(:,1), r_ode(:,2), r_ode(:,3)];
end
