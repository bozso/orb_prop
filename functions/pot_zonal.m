% Calculation of gravitational potential using zonal harmonics

function U = pot_zonal(r, n)
    global mu_si a_earth C

    delta = norm(r);
    rec_delta = 1.0 / delta;
    cos_theta = r(3) * rec_delta;
    
    a_per_r = a_earth * rec_delta;
    
    U = 1;
    
    for iii = 2:n
            Pn = legendre(iii, cos_theta);
            U = U + a_per_r^n * C(iii-1) * Pn(1) * sqrt((2*iii+1) * 0.5);
    end
    
    U = U * (mu_si * rec_delta);
end
