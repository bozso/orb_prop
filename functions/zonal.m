% Equation of motion with zonal harmonics.

function rdot = zonal(t, r)
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
