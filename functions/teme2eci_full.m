function t_poz_vel_eci = teme2eci_full(t_poz_vel_teme)
    global eod

    jday_UT1 = t_poz_vel_teme(:,1);

    [jd, jd_frac] = whole_and_frac(t_poz_vel_teme(:,1));
    
    dut1 = interp1(eod(:,1), eod(:,4), (jd + jd_frac) - 2400000.5);

    row = numel(jd);
    
    reci = zeros(row, 3);
    veci = zeros(row, 3);
    
    for iii = 1:row
        [year, mon, day, hr, min, sec] = invjday ( jd(iii), jd_frac(iii) );

        dat = iauDat(year, mon, day, jd_frac(iii));

        [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
         = convtime (year, mon, day, hr, min, sec, 0, dut1(iii), dat );
         
        [reci(iii,:), veci(iii,:), aeci] = teme2eci(t_poz_vel_teme(iii,2:4)', ...
                                                    t_poz_vel_teme(iii,5:end)', ...
                                                    [0 0 0]', ttt, 106, 2, 'a');
    end
    t_poz_vel_eci = [t_poz_vel_teme(:,1), reci, veci];
end
