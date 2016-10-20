function [t_poz_vel] = calc_theo(step, day, model, satrec)
    global whichconst day2sec sec2day eod zonal_end

    if (model == 'point_mass' || model == 'zonal_2' || model == 'zonal_2_4')
        % Getting initial conditions in TEME
        printf("Getting initial conditions from satrec...");
        [satrec, rteme ,vteme] = sgp4(satrec,  0.0);
        
        % Transforming into ECI
        
        [jd, jd_frac] = whole_and_frac(satrec.jdsatepoch);

        [year,mon,day,hr,min,sec] = invjday ( jd, jd_frac );

        dat = iauDat(year, mon, day, jd_frac);
        dut1 = interp1(eod(:,1), eod(:,4), satrec.jdsatepoch - 2400000.5);
        
        [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
             = convtime (year, mon, day, hr, min, sec, 0, dut1, dat );
        
        [reci, veci, aeci] = teme2eci(rteme', vteme', [0 0 0]', ttt, 106, 2, 'a');
        printf("DONE\n");
    end
	
    switch (model)
        case 'point_mass'
            disp('rk_propagate: Calculating with Earth as a point mass.');
        case 'zonal_2'
            printf("rk_propagate: J2 will be used, air ");
            printf("friction will be inored.\n");
            satrec.bstar = 0.0;
        case 'zonal_2_4'
            printf("rk_propagate: Zonal harmonics (J234) will be used, air ");
            printf("friction will be inored.\n");
            satrec.bstar = 0.0;
            sgp_teme = sgp4_propagate(step, day, model, satrec);
            sgp_eci = teme2eci(sgp_teme);
            sgp_eci = [sgp_eci(:,1), sgp_eci(:,2:end) * 1e3];
        otherwise
            disp('rk_propagate: Unrecognized model option!');
            return
    end

end
