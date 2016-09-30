% Calculates the orbital elements given the inital conditions: t_0, 
% location vector, velocity vector

function [a, ecc, incl, omega, sinw, cosw,  w, n, tau] = orbital_elements(t_0, rv_0)
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
            printf ("Warning: E_0 has a significant imaginary part: %f\n", E_0_im);
            printf ("Calculations, may be inaccurate.\n");
    else
            E_0 = real (E_0);
    end
    tau = t_0 - (1/n) * (E_0 - e * sin(E_0));
end
