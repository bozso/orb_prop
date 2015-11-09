function U = pot_zonal_j(r, theta)
	global mu_si a_earth J omega_earth

	delta = norm(r);
	cos_theta = r(3) / delta;
	
	a_per_r = a_earth / delta;
	
	U = 1;
	for n = 2:4
		Pn = legendre(n, cos_theta);
		U = U - a_per_r^n * J(n-1) * Pn(1) * sqrt(2*n+1);
	end
	
	U = U * (mu_si / delta);
end
