% Satellite trajectory calculation with ode78

function t_poz_vel = rk_propagate(step, day, model, jd0, rv0)
    % Processing parameters and initial conditions
    global sec2day day2sec omega

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

    % Propagation with ode78
    [t, r_ode] = ode78(f_handle, [0.0, day2sec * day], rv0, vopt);

    t_poz_vel = [(jd0 + t * sec2day), r_ode];
end
