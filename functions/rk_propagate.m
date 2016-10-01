% Satellite trajectory calculation with ode78

function [t_poz_vel] = rk_propagate(step, day, model, satrec)
    % Processing parameters and initial conditions
    global whichconst day2sec sec2day eod
	
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
    [satrec, rteme ,vteme] = sgp4(satrec,  0.0);
    
    % Transforming into ECI
    
    [jd, jd_frac] = whole_and_frac(satrec.jdsatepoch);

    [year,mon,day,hr,min,sec] = invjday ( jd, jd_frac );

    dat = iauDat(year, mon, day, jd_frac);
    dut1 = interp1(eod(:,1), eod(:,4), satrec.jdsatepoch - 2400000.5);
    
    [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
         = convtime (year, mon, day, hr, min, sec, 0, dut1, dat );
    
    [reci, veci, aeci] = teme2eci(rteme, vteme, [0 0 0], ttt, 106, 2, 'a');
    
    % Changing into SI units
    poz_vel_eci = [reci * 1e3, veci * 1e3];

    % Propagation with ode78
    [t, r_ode] = ode78(f_handle, [0.0, day2sec * day], poz_vel_eci', vopt);

    t_poz_vel = [(satrec.jdsatepoch + t .* sec2day), r_ode];
end
