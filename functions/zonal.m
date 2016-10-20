% Equation of motion with zonal harmonics.

function rdot = zonal(t, r, jd0)
    % Position vector
    rdot(1) = r(4);
    rdot(2) = r(5);
    rdot(3) = r(6);
    
	% Numerial derivation stepsize in meters.
    h = 1;
    
    [lod, xp, yp, jdut1, ttt] = get_transf_const(jd0, t);
    
    eqeterms = 2;
    ddpsi = ddeps = 0.0;

    % ---- find matrices
    [prec,psia,wa,ea,xa] = precess ( ttt, '80' );

    [deltapsi,trueeps,meaneps,omega,nut] = nutation(ttt,ddpsi,ddeps);

    [st,stdot] = sidereal(jdut1,deltapsi,meaneps,omega,lod,eqeterms );

    [pm] = polarm(xp,yp,ttt,'80');

    % ---- perform transformations
    thetasa= 7.29211514670698e-05 * (1.0  - lod/86400.0 );
    omegaearth = [0; 0; thetasa;];
    
    rpef  = st'*nut'*prec'*r(1:3);
    recef = pm'*rpef;

    g_ecef(1) = (pot_zonal(recef + [h, 0, 0]') - ... 
            pot_zonal(recef - [h, 0, 0]')) / (2 * h);

    g_ecef(2) = (pot_zonal(recef + [0, h, 0]') - ... 
            pot_zonal(recef - [0, h, 0]')) / (2 * h);

    g_ecef(3) = (pot_zonal(recef + [0, 0, h]') - ... 
            pot_zonal(recef - [0, 0, h]')) / (2 * h);

    g_pef = pm*g_ecef';
    g_eci = prec*nut*st*g_pef;
    
    rdot(4:6) = g_eci(1:3);
end
