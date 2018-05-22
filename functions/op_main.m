function out = op_main(fun_name, varargin)
    
    switch(fun_name)
        case 'sgp4_propagate'
            out = sgp4_propagate(varargin{:})
        case 'rk_propagate'
            out = rk_propagate(varargin{:})
        case 'plot_orbit'
            plot_orbit(varargin{:})
        otherwise
            error(['Unknown function ', fun_name]);
    end
end

function ellipse = calc_ellipse(t, tau, ecc, a, n)

    M = n * (t - tau);
    
    %for ii = 1 : numel(t)
            [E, sinv, cosv] = newtonm_scv (ecc, M);
    %end
    
    r = a * (1 - e * cos(E));
    ellipse(:,1) = r .* cosv;
    ellipse(:,2) = r .* sinv;
end

function out = calc_kepler(r0, v0, step, nsteps)

    r = zeros(nsteps, 3);
    v = zeros(nsteps, 3);

    r(1,:) = r0;
    v(1,:) = v0;

    for ii = 2:nsteps
        [r(ii,:), v(ii,:)] =  kepler(r(ii - 1,:) , v(ii - 1,:), step);
    end
    out.r = r;
    out.v = v;
end

function out = newtonm_scv ( ecc,m );
% Function from the SGP4/SDP4 library modified so it returns sin(v) and cos(v)
% instead of v.
% ------------------------------------------------------------------------------
%
%                           function newtonm
%
%  this function performs the newton rhapson iteration to find the
%    eccentric anomaly given the mean anomaly.  the true anomaly is also
%    calculated.
%
%  author        : david vallado                  719-573-2600    9 jun 2002
%
%  revisions
%                -
%
%  inputs          description                    range / units
%    ecc         - eccentricity                   0.0  to
%    m           - mean anomaly                   -2pi to 2pi rad
%
%  outputs       :
%    e0          - eccentric anomaly              0.0  to 2pi rad
%    nu          - true anomaly                   0.0  to 2pi rad
%
%  locals        :
%    e1          - eccentric anomaly, next value  rad
%    sinv        - sine of nu
%    cosv        - cosine of nu
%    ktr         - index
%    r1r         - cubic roots - 1 to 3
%    r1i         - imaginary component
%    r2r         -
%    r2i         -
%    r3r         -
%    r3i         -
%    s           - variables for parabolic solution
%    w           - variables for parabolic solution
%
%  coupling      :
%    cubic       - solves a cubic polynomial
%
%  references    :
%    vallado       2001, 72-75, alg 2, ex 2-1
%
% [e0,nu] = newtonm ( ecc,m );
% ------------------------------------------------------------------------------

        % -------------------------  implementation   -----------------
        numiter =    50;
        small   =     0.00000001;
        halfpi  = pi * 0.5;
        
        % -------------------------- hyperbolic  ----------------------
        if ( (ecc-1.0 ) > small )
           % -------------------  initial guess -----------------------
            if ( ecc < 1.6  )
                if ( ((m<0.0 ) & (m>-pi)) | (m>pi) )
                    e0= m - ecc;
                  else
                    e0= m + ecc;
                  end
              else
                if ( (ecc < 3.6 ) & (abs(m) > pi) )
                    e0= m - sign(m)*ecc;
                  else
                    e0= m/(ecc-1.0 );
                  end
              end
            ktr= 1;
            e1 = e0 + ( (m-ecc*sinh(e0)+e0) / (ecc*cosh(e0) - 1.0 ) );
            while ((abs(e1-e0)>small ) & ( ktr<=numiter ))
                e0= e1;
                e1= e0 + ( ( m - ecc*sinh(e0) + e0 ) / ( ecc*cosh(e0) - 1.0  ) );
                ktr = ktr + 1;
              end
            % ----------------  find true anomaly  --------------------
            sinv= -( sqrt( ecc*ecc-1.0  ) * sinh(e1) ) / ( 1.0  - ecc*cosh(e1) );
            cosv= ( cosh(e1) - ecc ) / ( 1.0  - ecc*cosh(e1) );
            nu  = atan2( sinv,cosv );
          else
            % --------------------- parabolic -------------------------
            if ( abs( ecc-1.0  ) < small )
%                c = [ 1.0/3.0; 0.0; 1.0; -m];
%                [r1r] = roots (c);
%                e0= r1r;
                 s = 0.5  * (halfpi - atan( 1.5 *m ) );
                 w = atan( tan( s )^(1.0 /3.0 ) );
                 e0= 2.0 *cot(2.0 *w);
                ktr= 1;
                nu = 2.0  * atan(e0);
              else
                % -------------------- elliptical ----------------------
                if ( ecc > small )
                    % -----------  initial guess -------------
                    if ( ((m < 0.0 ) & (m > -pi)) | (m > pi) )
                        e0= m - ecc;
                      else
                        e0= m + ecc;
                      end
                    ktr= 1;
                    e1 = e0 + ( m - e0 + ecc*sin(e0) ) / ( 1.0  - ecc*cos(e0) );
                    while (( abs(e1-e0) > small ) & ( ktr <= numiter ))
                        ktr = ktr + 1;
                        e0= e1;
                        e1= e0 + ( m - e0 + ecc*sin(e0) ) / ( 1.0  - ecc*cos(e0) );
                      end
                    % -------------  find true anomaly  ---------------
                    sinv= ( sqrt( 1.0 -ecc*ecc ) * sin(e1) ) / ( 1.0 -ecc*cos(e1) );
                    cosv= ( cos(e1)-ecc ) / ( 1.0  - ecc*cos(e1) );
                    %~ nu  = atan2( sinv,cosv );
                  else
                    % -------------------- circular -------------------
                    ktr= 0;
                    nu= m;
                    e0= m;
                  end
              end
          end
    
    out.e0   = e0
    out.sinv = sinv;
    out.cosv = cosv;
end

function Rot = calc_rotmatrix (sinw, cosw, incl, omega)
% Returns the Rot matrix that transforms coordinates from the orbital plane to the 
% reference plane.
    
    P_x = cosw * cos(omega) - sinw * sin(omega) * cos(incl);
    P_y = cosw * sin(omega) + sinw * cos(omega) * cos(incl);
    P_z = sinw * sin(incl);
    Q_x = -sinw * cos(omega) - cosw * sin(omega) * cos(incl);
    Q_y = -sinw * sin(omega) + cosw * cos(omega) * cos(incl);
    Q_z = cosw * sin(incl);

    Rot =  [P_x, Q_x; ...
            P_y, Q_y; ...
            P_z, Q_z]';
end

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
            fprintf('rk_propagate: Calculating with Earth as a point mass.\n');
        case 'zonal_2'
            printf("rk_propagate: J2 will be used, air ");
            printf("friction will be inored.\n");
            satrec.bstar = 0.0;
        case 'zonal_2_4'
            fprintf(["rk_propagate: Zonal harmonics (J234) will be used, air\n", ...
                     "friction will be inored.\n"]);
            satrec.bstar = 0.0;
            sgp_teme = sgp4_propagate(step, day, model, satrec);
            sgp_eci = teme2eci(sgp_teme);
            sgp_eci = [sgp_eci(:,1), sgp_eci(:,2:end) * 1e3];
        otherwise
            fprintf('rk_propagate: Unrecognized model option!\n');
            return
    end
end

function [delta] = compare(poz_1, poz_2)
% Compares two satellite trajectories
    if ( size(poz_1) ~= size(poz_2) )
        error('Size of the two data matrices is not equal!');
    end
    
    delta = poz_1 - poz_2;
    delta = sqrt(sum(delta.^2, 2));
end

function delta_sigma = deg_diff(r1, r2)

    row = rows(r1);

    delta_sigma = zeros(row, 1);

    [theta_1, phi_1, r_1] = cart2sph(r1);
    [theta_2, phi_2, r_2] = cart2sph(r2);

    for iii = 1:rows(r1)
        delta_sigma(iii) = polar_dist(theta_1(iii), theta_2(iii), ...
                                      phi_1(iii), phi_2(iii));
    end
end

function  deltat = iauDat(iy, im, id, fd)
%
%   - - - - - - -
%    i a u D a t
%   - - - - - - -
% 
%   For a given UTC date, calculate delta(AT) = TAI-UTC.
% 
%      :------------------------------------------:
%      :                                          :
%      :                 IMPORTANT                :
%      :                                          :
%      :  A new version of this function must be  :
%      :  produced whenever a new leap second is  :
%      :  announced.  There are four items to     :
%      :  change on each such occasion:           :
%      :                                          :
%      :  1) A new line must be added to the set  :
%      :     of statements that initialize the    :
%      :     array "changes".                     :
%      :                                          :
%      :  2) The parameter IYV must be set to     :
%      :     the current year.                    :
%      :                                          :
%      :  3) The "Latest leap second" comment     :
%      :     below must be set to the new leap    :
%      :     second date.                         :
%      :                                          :
%      :  4) The "This revision" comment, later,  :
%      :     must be set to the current date.     :
%      :                                          :
%      :  Change (2) must also be carried out     :
%      :  whenever the function is re-issued,     :
%      :  even if no leap seconds have been       :
%      :  added.                                  :
%      :                                          :
%      :  Latest leap second:  2012 June 30       :
%      :                                          :
%      :__________________________________________:
% 
%   This function is part of the International Astronomical Union's
%   SOFA (Standards Of Fundamental Astronomy) software collection.
% 
%   Status:  support function.
% 
%   Given:
%      iy     int      UTC:  year (Notes 1 and 2)
%      im     int            month (Note 2)
%      id     int            day (Notes 2 and 3)
%      fd     double         fraction of day (Note 4)
% 
%   Returned:
%      deltat double   TAI minus UTC, seconds
% 
%   Returned (function value):
%             int      status (Note 5):
%                        1 = dubious year (Note 1)
%                        0 = OK
%                       -1 = bad year
%                       -2 = bad month
%                       -3 = bad day (Note 3)
%                       -4 = bad fraction (Note 4)
% 
%   Notes:
% 
%   1) UTC began at 1960 January 1.0 (JD 2436934.5) and it is improper
%      to call the function with an earlier date.  If this is attempted,
%      zero is returned together with a warning status.
% 
%      Because leap seconds cannot, in principle, be predicted in
%      advance, a reliable check for dates beyond the valid range is
%      impossible.  To guard against gross errors, a year five or more
%      after the release year of the present function (see parameter
%      IYV) is considered dubious.  In this case a warning status is
%      returned but the result is computed in the normal way.
% 
%      For both too-early and too-late years, the warning status is
%      j=+1.  This is distinct from the error status j=-1, which
%      signifies a year so early that JD could not be computed.
% 
%   2) If the specified date is for a day which ends with a leap second,
%      the UTC-TAI value returned is for the period leading up to the
%      leap second.  If the date is for a day which begins as a leap
%      second ends, the UTC-TAI returned is for the period following the
%      leap second.
% 
%   3) The day number must be in the normal calendar range, for example
%      1 through 30 for April.  The "almanac" convention of allowing
%      such dates as January 0 and December 32 is not supported in this
%      function, in order to avoid confusion near leap seconds.
% 
%   4) The fraction of day is used only for dates before the
%      introduction of leap seconds, the first of which occurred at the
%      end of 1971.  It is tested for validity (0 to 1 is the valid
%      range) even if not used;  if invalid, zero is used and status
%      j=-4 is returned.  For many applications, setting fd to zero is
%      acceptable;  the resulting error is always less than 3 ms (and
%      occurs only pre-1972).
% 
%   5) The status value returned in the case where there are multiple
%      errors refers to the first error detected.  For example, if the
%      month and day are 13 and 32 respectively, j=-2 (bad month)
%      will be returned.
% 
%   6) In cases where a valid result is not available, zero is returned.
% 
%   References:
% 
%   1) For dates from 1961 January 1 onwards, the expressions from the
%      file ftp://maia.usno.navy.mil/ser7/tai-utc.dat are used.
% 
%   2) The 5ms timestep at 1961 January 1 is taken from 2.58.1 (p87) of
%      the 1992 Explanatory Supplement.
% 
%   Called:
%      iauCal2jd    Gregorian calendar to Julian Day number
% 
%   This revision:  2012 June 14
% 
%   SOFA release 2012-03-01
% 
%   Copyright (C) 2012 IAU SOFA Board.  See notes at end.
%

    % Release year for this version of iauDat
    IYV = 2012;
    
    % Reference dates (MJD) and drift rates (s/day), pre leap seconds
    drift = [
       37300.0, 0.0012960 ;
       37300.0, 0.0012960 ;
       37300.0, 0.0012960 ;
       37665.0, 0.0011232 ;
       37665.0, 0.0011232 ;
       38761.0, 0.0012960 ;
       38761.0, 0.0012960 ;
       38761.0, 0.0012960 ;
       38761.0, 0.0012960 ;
       38761.0, 0.0012960 ;
       38761.0, 0.0012960 ;
       38761.0, 0.0012960 ;
       39126.0, 0.0025920 ;
       39126.0, 0.0025920
    ];
    
    %  Number of Delta(AT) expressions before leap seconds were introduced 
    [NERA1 temp] = size (drift);
    
     % iyear month  Delta(AT)
    changes = [
       1960,  1,  1.4178180 ;
       1961,  1,  1.4228180 ;
       1961,  8,  1.3728180 ;
       1962,  1,  1.8458580 ;
       1963, 11,  1.9458580 ;
       1964,  1,  3.2401300 ;
       1964,  4,  3.3401300 ;
       1964,  9,  3.4401300 ;
       1965,  1,  3.5401300 ;
       1965,  3,  3.6401300 ;
       1965,  7,  3.7401300 ;
       1965,  9,  3.8401300 ;
       1966,  1,  4.3131700 ;
       1968,  2,  4.2131700 ;
       1972,  1, 10.0       ;
       1972,  7, 11.0       ;
       1973,  1, 12.0       ;
       1974,  1, 13.0       ;
       1975,  1, 14.0       ;
       1976,  1, 15.0       ;
       1977,  1, 16.0       ;
       1978,  1, 17.0       ;
       1979,  1, 18.0       ;
       1980,  1, 19.0       ;
       1981,  7, 20.0       ;
       1982,  7, 21.0       ;
       1983,  7, 22.0       ;
       1985,  7, 23.0       ;
       1988,  1, 24.0       ;
       1990,  1, 25.0       ;
       1991,  1, 26.0       ;
       1992,  7, 27.0       ;
       1993,  7, 28.0       ;
       1994,  7, 29.0       ;
       1996,  1, 30.0       ;
       1997,  7, 31.0       ;
       1999,  1, 32.0       ;
       2006,  1, 33.0       ;
       2009,  1, 34.0       ;
       2012,  7, 35.0
    ];
    
    %  Number of Delta(AT) changes 
    NDAT = size(changes);
    
    %  Initialize the result to zero. 
    deltat = 0.0;
    da = 0.0;
    
    %  If invalid fraction of a day, set error status and give up. 
    if (fd < 0.0 || fd > 1.0) 
        error('invalid fraction of a day')
    end
    
    %  Convert the date into an MJD. 
    [djm0 djm] = iauCal2jd(iy, im, id);
    
    %  If pre-UTC year, set warning status and give up. 
    if (iy < changes(1,1))
        warning('pre-UTC year')
        exit
    end
    
    %  If suspiciously late year, set warning status but proceed. 
    if (iy > IYV + 5) 
        warning('suspiciously late year')
    end
    
    %  Combine year and month to form a date-ordered integer... 
    m = 12*iy + im;
    
    %  ...and use it to find the preceding table entry.
    for i = NDAT:-1:1
        if (m >= (12 * changes(i,1) + changes(i,2))) 
            break
        end
    end
    
    %  Get the Delta(AT)
    da = changes(i,3);
    
    %  If pre-1972, adjust for drift
    if (i < NERA1)
        da = da + (djm + fd - drift(i,1)) * drift(i,2);
    end
    
    %  Return the Delta(AT) value
    deltat = da;
end

function [djm0 djm] = iauCal2jd(Year, Month, Day, Hour, Min, Sec)
%--------------------------------------------------------------------------
% Purpose:
%
%   Modified Julian Date from calendar date and time
%
% Input(s):
%
%   Year      
%   Month
%   Day
%   Hour      
%   Min
%   Sec
% 
% output(s):  
%    
%    djm0  
%    djm     Modified Julian Date
% 
% Reference: O. Montenbruck, E. Gill (2000), "Satellite Orbits"
% 
%--------------------------------------------------------------------------

    if (nargin == 3)
        Hour = 0;
        Min = 0;
        Sec = 0;
    end
    
    djm0 = 2400000.5;
    
    if ( Month<=2 ) 
        Month = Month+12;
        Year = Year-1;
    end
    
    if ( (10000*Year+100*Month+Day) <= 15821004 )
        b = fix(-2 + ((Year+4716)/4) - 1179);     % Julian calendar
    else
        b = fix((Year/400)-(Year/100)+(Year/4));  % Gregorian calendar
    end
    
    MjdMidnight = 365*Year - 679004 + b + fix(30.6001*(Month+1)) + Day;
    FracOfDay   = (Hour+Min/60.0+Sec/3600.0) / 24.0;
    
    djm = MjdMidnight + FracOfDay;
end

function xlabel = get_xlabel(j0)

    [jd, jdfrac] = whole_and_frac(j0);

    [year, mon, day, hr, min, sec] = invjday( jd, jdfrac );

    xlabel = ['Days from ', num2str(year), '.', double_digit(mon), '.', ...
              double_digit(day), ' ', double_digit(hr), ':', ...
              double_digit(min), ':', double_digit(sec)];

end

function out = get_transf_const(jd0, t)
    global eod
    
    [jd, jd_frac] = whole_and_frac(jd0 + t / 86400);
    
    jdut1 = jd + jd_frac;
    
    dut1 = interp1(eod(:,1), eod(:,4), (jd + jd_frac) - 2400000.5);
    lod = interp1(eod(:,1), eod(:,4), (jd + jd_frac) - 2400000.5);
    xp = interp1(eod(:,1), eod(:,2), (jd + jd_frac) - 2400000.5) ./ 206265;
    yp = interp1(eod(:,1), eod(:,3), (jd + jd_frac) - 2400000.5) ./ 206265;

    [year, mon, day, hr, min, sec] = invjday ( jd, jd_frac );
    
    dat = iauDat(year, mon, day, jd_frac);
    
    [ut1, tut1, jdut1,jdut1frac, utc, tai, tt, ttt, jdtt,jdttfrac, tdb, ttdb, ...
     jdtdb,jdtdbfrac, tcg, jdtcg,jdtcgfrac, tcb, jdtcb,jdtcbfrac ] ...
     = convtime (year, mon, day, hr, min, sec, 0, dut1, dat );
    
    out.lod     = lod;
    out.xp      = xp;
    out.yp      = yp;
    out.jd      = jd
    out.jd_frac = jd_frac;
    out.ttt     = ttt;
end

function rv0 = get_init(satrec)

    % Getting initial conditions in TEME
    fprintf("Getting initial conditions from satrec...\n");
    [satrec, rteme ,vteme] = sgp4(satrec,  0.0);

    % Transforming into ECEF
    t_poz_vel_0 = teme2ecef_full([satrec.jdsatepoch, rteme ,vteme]);

    % Changing into SI units
    rv0 = [t_poz_vel_0(2:4) * 1e3, t_poz_vel_0(5:end) * 1e3];
end

function digit = double_digit(num)
    if(num < 10)
        digit = ['0', num2str(num)];
    else
        digit = num2str(num);
    end
end

function [a, ecc, incl, omega, sinw, cosw,  w, n, tau] = orbital_elements(t_0, rv_0)
% Calculates the orbital elements given the inital conditions: t_0, 
% location vector, velocity vector
    global mu_si
    r_0 = [rv_0(1), rv_0(2), rv_0(3)];
    v_0 = [rv_0(4), rv_0(5), rv_0(6)];
    
    c = cross(r_0, v_0);
    r = norm(r_0);
    per_r0 = 1 / r;

    lambda = ( dot(v_0, v_0) - (mu_si * per_r0) ) * r_0 + - dot(r_0, v_0) * v_0;
    
    dot (c, lambda)
    
    incl = atan2( norm ( [c(1) c(2)] ), c(3) );
    
    omega = atan2(-c(1), c(2)) + pi;
    
    norm_lambda = norm(lambda);
    
    sinw = ( -lambda(1) * sin(omega) + lambda(2) * cos(omega) ) ...
            / ( norm_lambda * cos(incl) );
    cosw = ( lambda(1) * cos(omega) + lambda(2) * sin(omega) ) ...
            / norm_lambda;
    
    w = atan2 ( -lambda(1) * sin(omega) + lambda(2) * cos(omega), ...
                lambda(1) * cos(omega) + lambda(2) * sin(omega) );
    
    h = (v_0 * v_0') / 2 - (mu_si * per_r0);
    
    a = - mu_si / (2*h);
    per_a = 1 / a;
    n = sqrt( mu_si * (per_a^3) );
    ecc = sqrt( 1 + 2 * h * ( dot(c, c) / (mu_si * mu_si) ) );
    E_0 = acos( ( 1 - r * per_a) / e );
    E_0_im = imag(E_0);
    
    if ( E_0_im > 1e-3)
        warning(["E_0 has a significant imaginary part: %f\n",
                 "Calculations, may be inaccurate.\n"], E_0_im);
    else
        E_0 = real (E_0);
    end
    tau = t_0 - (1/n) * (E_0 - e * sin(E_0));
end

function rdot = zonal(t, r)
% Equation of motion with zonal harmonics.
    global omega
    
    % Position vector   
    rdot(1:3) = r(4:6);
    
	% Numerial derivation stepsize in meters.
    h = 1;
    
    g_ecef(1) = (pot_zonal(r(1:3) + [h, 0, 0]') - ... 
            pot_zonal(r(1:3) - [h, 0, 0]')) / (2 * h);

    g_ecef(2) = (pot_zonal(r(1:3) + [0, h, 0]') - ... 
            pot_zonal(r(1:3) - [0, h, 0]')) / (2 * h);

    g_ecef(3) = (pot_zonal(r(1:3) + [0, 0, h]') - ... 
            pot_zonal(r(1:3) - [0, 0, h]')) / (2 * h);
        
    rdot(4:6) = g_ecef' - 2 * cross(omega, r(4:6)) - cross(omega, ... 
                cross(omega, r(1:3)));
end

function [whole, frac] = whole_and_frac(num)
	whole = floor(num);
	frac = num - whole;
end

function [t_poz_vel] = sgp4_propagate(step, day, model, satrec)
% Satellite trajectory calculation with SGP

    global whichconst day2sec sec2day
	
	time = 0:step:day*day2sec;
	
	t_poz_vel = zeros(numel(time), 7);
	
    switch (model)
        case 'point_mass'
            error(["sgp4_propagate: SGP4 can not handle the Earth as ", ...
                   "a point mass.\n"]);
        case 'zonal'
            fprintf(["sgp4_propagate: Zonal harmonics will be used, air ", ...
                     "friction will be inored.\n"]);
            satrec.bstar = 0.0;
        otherwise
            error('Unrecognized model option!');
    end
    
    global tumin mu radiusearthkm xke j2 j3 j4 j3oj2 opsmode

	sec2min = 1 / 60;

    for ii = 1:numel(time)
        [satrec, ro ,vo] = sgp4 (satrec,  time(ii) * sec2min);

        if (satrec.error > 0)
            fprintf('# *** error: t:= %f *** code = %3i\n', time(ii), ...
                    satrec.error);
        end  
                
        if (satrec.error == 0)
            t_poz_vel(iii,:) = [satrec.jdsatepoch + time(ii) * sec2day,...
                                ro, vo];
%            [satrec.jdsatepoch + time(iii) * sec2day, ro, vo]
%            t_poz_vel(iii,:)
        else
            error("sgp4_propagate: Error: satrec.error is non zero!");
        end
    end %// while propagating the orbit   
end	

function t_poz_vel = rk_propagate(step, day, model, jd0, rv0)
% Satellite trajectory calculation with ode78

    % Processing parameters and initial conditions
    global sec2day day2sec omega

    switch (model)
        case 'point_mass'
            fprintf('rk_propagate: Calculating with Earth as a point mass.\n');
            f_handle = @point_mass;
        case 'zonal'
            f_handle = @zonal;
            fprintf(["rk_propagate: Zonal harmonics (J234) will be used, air ", ...
                     "friction will be inored.\n"]);
        otherwise
            error('Unrecognized model option!');
    end

    vopt = odeset('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', ...
                  'InitialStep', step, 'MaxStep', step);

    % Propagation with ode78
    [t, r_ode] = ode78(f_handle, [0.0, day2sec * day], rv0, vopt);

    t_poz_vel = [(jd0 + t * sec2day), r_ode];
end

function [step day model satrec] = process_infile(file_loc)
% Processes inputfile

    global whichconst jday_low jday_high
    
    % File that contains propagation parameters and initial conditions
    infile = fopen(file_loc, 'r');
    
    % Processing parameters and initial conditions
    
    str = strsplit(fgetl(infile));
    
    if (numel(str) > 3)
        error('Too many options given in the first line!')
    end
    
    step = str2num(str{1});
    day = str2num(str{2});
    model = str{3};
    
    % Processing oribtal element set
    fprintf('Reading satellite data...')
    longstr1 = fgets(infile, 130);
    while ( (longstr1(1) == '#') && (feof(infile) == 0) )
            longstr1 = fgets(infile, 130);
    end
    
    if (feof(infile) == 0)
            longstr2 = fgets(infile, 130);
    end
            
#    satrec = process(whichconst, longstr1, longstr2);
#    [satrec, startmfe, stopmfe, deltamin] = twoline2rv(whichconst, longstr1, ...
#          longstr2, 'c','e');
    [startmfe, stopmfe, deltamin, satrec] = twoline2rv(longstr1, longstr2, ...
          'c', 'm', 'i', whichconst);
    fclose(infile);
    fprintf('DONE\n');

    if (satrec.jdsatepoch > jday_high || satrec.jdsatepoch < jday_low)
        warning(['Can only compare sgp and Runge-Kutta if ephemeris date ',...
                 'is between 1992.01.01 2016.07.01']);
	end
end

function satrec = process(whichconst, longstr1, longstr2)
% Process tle

    global tumin radiusearthkm xke j2 j3 j4 j3oj2  

    deg2rad  =   pi / 180.0;         %  0.01745329251994330;  % [deg/rad]
    xpdotp   =  1440.0 / (2.0*pi);   % 229.1831180523293;  % [rev/day]/[rad/min]  

    revnum = 0; 
    elnum  = 0;
    year   = 0; 
    satrec.error = 0;

%     // set the implied decimal points since doing a formated read
%     // fixes for bad input data values (missing, ...)
    for (j = 11:16)
        if (longstr1(j) == ' ')
            longstr1(j) = '_';
        end
    end

    if (longstr1(45) ~= ' ')
        longstr1(44) = longstr1(45);
    end
    longstr1(45) = '.';
     
    if (longstr1(8) == ' ')
        longstr1(8) = 'U';
    end

    if (longstr1(10) == ' ')
        longstr1(10) = '.';
    end

    for (j = 46:50)
        if (longstr1(j) == ' ')
            longstr1(j) = '0';
        end
    end
    if (longstr1(52) == ' ')
        longstr1(52) = '0';
    end
    if (longstr1(54) ~= ' ')
        longstr1(53) = longstr1(54);
    end
    longstr1(54) = '.';

    longstr2(26) = '.';
     
    for (j = 27:33)
        if (longstr2(j) == ' ')
            longstr2(j) = '0';
        end
    end
     
    if (longstr1(63) == ' ')
        longstr1(63) = '0';
    end

    if ((length(longstr1) < 68) || (longstr1(68) == ' '))
        longstr1(68) = '0';
    end

    % parse first line
    carnumb = str2num(longstr1(1));
    satrec.satnum = str2num(longstr1(3:7));
    classification = longstr1(8);
    intldesg = longstr1(10:17);
    satrec.epochyr = str2num(longstr1(19:20));
    satrec.epochdays = str2num(longstr1(21:32));
    satrec.ndot = str2num(longstr1(34:43));
    satrec.nddot = str2num(longstr1(44:50));
    nexp = str2num(longstr1(51:52));
    satrec.bstar = str2num(longstr1(53:59));
    ibexp = str2num(longstr1(60:61));
    numb = str2num(longstr1(63));
    elnum = str2num(longstr1(65:68));
 
    % parse second line
    cardnumb = str2num(longstr2(1));
    satrec.satnum = str2num(longstr2(3:7));
    satrec.inclo = str2num(longstr2(8:16));
    satrec.nodeo = str2num(longstr2(17:25));
    satrec.ecco = str2num(longstr2(26:33));
    satrec.argpo = str2num(longstr2(34:42));
    satrec.mo = str2num(longstr2(43:51));
    satrec.no = str2num(longstr2(52:63));
    revnum = str2num(longstr2(64:68));

%     // ---- find no, ndot, nddot ----
    satrec.no   = satrec.no / xpdotp; %//* rad/min
    satrec.nddot= satrec.nddot * 10.0^nexp;
    satrec.bstar= satrec.bstar * 10.0^ibexp;

%     // ---- convert to sgp4 units ----
    satrec.a    = (satrec.no*tumin)^(-2/3);                % [er]
    satrec.ndot = satrec.ndot  / (xpdotp*1440.0);          % [rad/min^2]
    satrec.nddot= satrec.nddot / (xpdotp*1440.0*1440);     % [rad/min^3]

%     // ---- find standard orbital elements ----
    satrec.inclo = satrec.inclo  * deg2rad;
    satrec.nodeo = satrec.nodeo * deg2rad;
    satrec.argpo = satrec.argpo  * deg2rad;
    satrec.mo    = satrec.mo     *deg2rad;

    satrec.alta = satrec.a*(1.0 + satrec.ecco) - 1.0;
    satrec.altp = satrec.a*(1.0 - satrec.ecco) - 1.0;

%     // ----------------------------------------------------------------
%     // find sgp4epoch time of element set
%     // remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
%     // and minutes from the epoch (time)
%     // --------------------------------------------------------------

%     // ------------- temp fix for years from 1957-2056 ----------------
%     // ------ correct fix will occur when year is 4-digit in 2le ------
     if (satrec.epochyr < 57)
         year= satrec.epochyr + 2000;
       else
         year= satrec.epochyr + 1900;
     end;

     [mon,day,hr,minute,sec] = days2mdh ( year,satrec.epochdays );
     satrec.jdsatepoch = jday( year,mon,day,hr,minute,sec );

     
%     // ------------- initialize the orbit at sgp4epoch --------------
     sgp4epoch = satrec.jdsatepoch - 2433281.5; % days since 0 Jan 1950
     [satrec] = sgp4init(whichconst, satrec, satrec.bstar, satrec.ecco, sgp4epoch, ...
         satrec.argpo, satrec.inclo, satrec.mo, satrec.no, satrec.nodeo);
end

function U = pot_zonal(r)
% Calculation of gravitational potential using zonal harmonics
    
    global mu_si a_earth C

    delta = norm(r);
    rec_delta = 1.0 / delta;
    cos_theta = r(3) * rec_delta;
    
    a_per_r = a_earth * rec_delta;
    
    U = 1;
    
    for n = 2:4
            Pn = legendre(n, cos_theta);
            U = U + a_per_r^n * C(n-1) * Pn(1) * sqrt((2*n+1) * 0.5);
    end
    
    U = U * mu_si * rec_delta;
end

function delta_sigma = polar_dist(theta_1, theta_2, phi_1, phi_2)
   delta_theta = abs(theta_1 - theta_2); 
   delta_phi= abs(phi_1 - phi_2); 

   cos_phi_1 = cos(phi_1);
   cos_phi_2 = cos(phi_2);
   sin_phi_1 = sin(phi_1);
   sin_phi_2 = sin(phi_2);

   cos_delta_theta = cos(delta_theta);
   
   a = (cos_phi_2 * sin(delta_theta))^2;
   b = (cos_phi_1 * sin_phi_2 - sin_phi_1 * cos_phi_2 * cos_delta_theta)^2;
   c = sin_phi_1 * sin_phi_2 + cos_phi_1 * cos_phi_2 * cos_delta_theta;

   delta_sigma = atan2( sqrt(a + b), c);
end

function rdot = point_mass(t, r)
% Equations of motion, assuming that the Earth is a point mass

    global mu_si
    delta = norm([r(1), r(2), r(3)]);
    mu_per_rcubed = mu_si / delta^3;
    rdot(1) = r(4);
    rdot(2) = r(5);
    rdot(3) = r(6);
    rdot(4) = - mu_per_rcubed * r(1);
    rdot(5) = - mu_per_rcubed * r(2);
    rdot(6) = - mu_per_rcubed * r(3);
end

function [] = plot_orbit(xyz, range, unit)
% Satellite trajectory plot

    plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3));

    if (numel(range) > 1)
        axis(range);
    end
    view(3);
    xlabel(['x [' unit ']']);
    ylabel(['y [' unit ']']);
    zlabel(['z [' unit ']']);
end
