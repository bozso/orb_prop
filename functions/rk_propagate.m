% Satellite trajectory calculation with ode78

function [t_poz_vel] = rk_propagate(step, day, model, satrec)
    % Processing parameters and initial conditions
    global whichconst day2sec
	
    switch (model)
        case 'point_mass'
            disp('rk_propagate: Calculating with Earth as a point mass.');
            f_handle = @point_mass;
        case 'zonal'
            f_handle = @zonal;
            printf("sgp4_propagate: Zonal harmonics will be used, air ");
            printf("friction will be inored.\n");
            satrec.bstar = 0.0;
        otherwise
            disp('rk_propagate: Unrecognized model option!');
            return
    end

    vopt = odeset('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', ...
					'InitialStep', 100, 'MaxStep', 100);

    % Getting initial conditions
    
    % Transformation!!!!!
    
    [satrec, ro ,vo] = sgp4(satrec,  0.0);

    % Changing into SI units
    rv_0 = [ro, vo] * 1e3;

    % Propagation with ode78
    [t, r_ode] = ode78(f_handle, [0.0, day2sec * day], rv_0, vopt);

    t_poz_vel = [(satrec.jdsatepoch + t .* sec2day)', r_ode];
end
