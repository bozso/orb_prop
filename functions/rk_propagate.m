% Satellite trajectory calculation with ode78

function [timexyzv] = rk_propagate(infile_path, outname)
    % Processing parameters and initial conditions
    global whichconst
	
	[step day model satrec] = process_infile(infile_path, whichconst);
	
    switch (model)
            case 'point_mass'
				disp('rk_propagate: Calculating with Earth as a point mass.');
				f_handle = @point_mass;
            case 'zonal'
				f_handle = @zonal;
				disp('rk_propagate: Zonal harmonics will be used, air friction ...
						will be inored')
				satrec.bstar = 0.0;
            otherwise
				disp('rk_propagate: Unrecognized model option!');
				return
    end

    vopt = odeset('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', ...
					'InitialStep', 100, 'MaxStep', 100);

    % Getting initial conditions
    [satrec, ro ,vo] = sgp4(satrec,  0.0);

    % Changing into SI units
    rv_0 = [ro, vo] * 1e3;

    % plotname is satellite number
    plotname = ['output/', outname];

    % Propagation with ode78
    [t, r_ode] = ode78(f_handle, [0.0, 86400 * day], rv_0, vopt);

    timexyzv = [t, r_ode];

    save('-ascii', ['output/', outname, '_rk.dat'], 'timexyzv');
end
