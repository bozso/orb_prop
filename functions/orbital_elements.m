%~ Calculates the orbital elements given the inital conditions: t_0, location vector, velocity vector

function [a, e, i, omega, sinw, cosw,  n, tau] = orbital_elements(t_0, rv_0)
	global mu_si
	r_0 = [rv_0(1), rv_0(2), rv_0(3)];
	v_0 = [rv_0(4), rv_0(5), rv_0(6)];
	
	c = cross(r_0, v_0);
	r = norm(r_0);
	
	lambda = - (mu_si / r)* r_0 + cross (v_0, c);
	
	i = - atan2(sqrt(c(1) * c(1) + c(2) * c(2)), c(3));
	
	if (-c(1) == 0 && c(2) == 0)
		omega = 0.0;
	else
		omega = atan2(-c(1), c(2));
	end
	
	sinw = ( -lambda(1) * sin(omega) + lambda(2) * cos(omega) ) / ( norm(lambda) * cos(i) );	
	cosw = ( lambda(2) * sin(omega) + lambda(1) * cos(omega) ) / ( norm(lambda) );

	h = (v_0 * v_0') / 2 - (mu_si / r);
	
	a = -mu_si / (2*h);
	n = sqrt(mu_si / (a*a*a));
	e = sqrt( 1 + 2 * h * (c*c' / (mu_si * mu_si)) );
	E_0 = acos( (1 - (r / a)) / e );
	
	E_O_im = imag(E_0);
	
	if ( E_O_im > 1e-3)
		printf ('Warning: E_0 has a significant imaginary part: %f\n', E_0_im);
		printf ('Calculations, may be inaccurate.\n');
	else
		E_0 = real (E_0);
	end
	tau = t_0 - (1/n) * (E_0 - e * sin(E_0));
end
