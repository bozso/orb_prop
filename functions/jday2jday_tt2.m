function [] = jday2jday_tt2.m()
    [jd, jdfrac] =  jday( year,mon,day,0,0,0.0 );
    mjd  = jd+jdfrac - 2400000.5;
    mfme = hr*60.0 + min + sec/60.0;

    % ------------------ start if ( ut1 is known ------------------
    localhr= timezone + hr;

    utc = hms2sec( localhr,min,sec );

    ut1= utc + dut1;
    [hrtemp,mintemp,sectemp] = sec2hms(  ut1 );
    [jdut1, jdut1frac] = jday( year,mon,day, hrtemp, mintemp, sectemp );
    tut1= (jdut1+jdut1frac - 2451545.0  )/ 36525.0;

    tai= utc + dat;

    tt= tai + 32.184;   % sec
    [hrtemp,mintemp,sectemp] = sec2hms( tt );
    [jdtt, jdttfrac] = jday( year,mon,day, hrtemp, mintemp, sectemp);
    ttt= (jdtt - 2451545.0  )/ 36525.0;
end
