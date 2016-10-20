% Satellite trajectory calculation with ode78

function t_poz_vel = rk_propagate(step, day, model, satrec)
    % Processing parameters and initial conditions
    global whichconst day2sec sec2day eod zonal_end jd0
	
    jd0 = satrec.jdsatepoch;
    
    switch (model)
        case 'point_mass'
            disp('rk_propagate: Calculating with Earth as a point mass.');
            f_handle = @point_mass;
        case 'zonal'
            f_handle = @zonal;
            printf("rk_propagate: Zonal harmonics (J234) will be used, air ");
            printf("friction will be inored.\n");
        otherwise
            disp('rk_propagate: Unrecognized model option!');
            return
    end

    vopt = odeset('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', ...
					'InitialStep', step, 'MaxStep', step);

    % Getting initial conditions in TEME
    printf("Getting initial conditions from satrec...");
    [satrec, rteme ,vteme] = sgp4(satrec,  0.0);
    
    % Transforming into ECI
    
    t_poz_vel_0 = teme2eci_full([jd0, rteme ,vteme]);
        
    % Changing into SI units
    poz_vel_eci = [t_poz_vel_0(2:4) * 1e3, t_poz_vel_0(5:end) * 1e3];
    
    printf("DONE\n");
    
    % Propagation with ode78
    printf("Propagating with ode87...");
    [t, r_ode] = ode78(f_handle, [0.0, day2sec * day], poz_vel_eci', vopt, ...
                        jd0);
    printf("DONE\n");
    
    t_poz_vel = [(jd0 + t .* sec2day), r_ode];
end
