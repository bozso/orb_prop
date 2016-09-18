% Satellite trajectory calculation with ode78

function [t_poz_vel] = rk_propagate(step, day, model, satrec)
    % Processing parameters and initial conditions
    global whichconst day2sec sec2day
	
    switch (model)
        case 'point_mass'
            disp('rk_propagate: Calculating with Earth as a point mass.');
            f_handle = @point_mass;
        case 'zonal'
            f_handle = @zonal;
            printf("rk_propagate: Zonal harmonics will be used, air ");
            printf("friction will be inored.\n");
            satrec.bstar = 0.0;
        otherwise
            disp('rk_propagate: Unrecognized model option!');
            return
    end

    vopt = odeset('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', ...
					'InitialStep', 100, 'MaxStep', 100);

    % Getting initial conditions in TEME
    [satrec, r_teme ,v_teme] = sgp4(satrec,  0.0);
    
    % Transforming into ECI
    poz_vel_eci = teme2eci([satrec.jdsatepoch, r_teme, v_teme]);
    
    % Changing into SI units
    poz_vel_eci = poz_vel_eci(:,2:end) * 1e3;

    % Propagation with ode78
    [t, r_ode] = ode78(f_handle, [0.0, day2sec * day], poz_vel_eci', vopt);

    t_poz_vel = [(satrec.jdsatepoch + t .* sec2day), r_ode];
end
