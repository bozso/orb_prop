function [lod, xp, yp, jd, jd_frac, ttt]  = get_transf_const(jd0, t)
    global eod
    
    [jd, jd_frac] = whole_and_frac(jd0 + t / 86400);
    
    jdut1 = jd + jd_frac;
    
    dut1 = interp1(eod(:,1), eod(:,4), (jd + jd_frac) - 2400000.5);
    lod = interp1(eod(:,1), eod(:,4), (jd + jd_frac) - 2400000.5);
    xp = interp1(eod(:,1), eod(:,2), (jd + jd_frac) - 2400000.5) ./ 206265;
    yp = interp1(eod(:,1), eod(:,3), (jd + jd_frac) - 2400000.5) ./ 206265;

    [year, mon, day, hr, min, sec] = invjday ( jd, jd_frac );
    
    dat = iauDat(year, mon, day, jd_frac);
    
    [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
     = convtime (year, mon, day, hr, min, sec, 0, dut1, dat );
end
