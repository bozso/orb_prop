function U = pot_zonal(r)
    global mu_si a_earth C

    delta = norm(r);
    rec_delta = 1.0 / delta;
    cos_theta = r(3) * rec_delta;
    
    a_per_r = a_earth * rec_delta;
    
    U = 1;
    for n = 2:4
            Pn = legendre(n, cos_theta);
            U = U + a_per_r^n * C(n-1) * Pn(1) * sqrt(2*n+1);
    end
    
    U = U * (mu_si * rec_delta);
end
