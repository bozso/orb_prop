% transform

function t_poz_vel_ecef = eci2ecef_full(t_poz_vel_eci)
    global eod

    [jd, jd_frac] = whole_and_frac(t_poz_vel_eci(:,1));
    
    dut1 = interp1(eod(:,1), eod(:,4), (jd + jd_frac) - 2400000.5);
    lod = interp1(eod(:,1), eod(:,4), (jd + jd_frac) - 2400000.5);
    xp = interp1(eod(:,1), eod(:,2), (jd + jd_frac) - 2400000.5) ./ 206265;
    yp = interp1(eod(:,1), eod(:,3), (jd + jd_frac) - 2400000.5) ./ 206265;

    row = numel(jd);
    
    recef = zeros(row, 3);
    vecef = zeros(row, 3);
    
    for iii = 1:row
        [year, mon, day, hr, min, sec] = invjday ( jd(iii), jd_frac(iii) );

        dat = iauDat(year, mon, day, jd_frac(iii));

        [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
         = convtime (year, mon, day, hr, min, sec, 0, dut1(iii), dat );
        
        [recef(iii,:), vecef(iii,:), aeci] = eci2ecef(t_poz_vel_eci(iii,2:4)', ...
                                    t_poz_vel_eci(iii,5:end)', ...
                                    [0 0 0]', ttt, jd(iii) + jd_frac(iii),...
                                    lod(iii), xp(iii), yp(iii), 2, 0.0, 0.0);
    end
    t_poz_vel_ecef = [t_poz_vel_eci(:,1), recef, vecef];
end
