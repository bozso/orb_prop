function rdot = Gravi(t, r)
	GRAVITY_CONST = 6.67384e-11;
	EARTH_MASS = 5.97237e24;
	delta = sqrt(r(1)*r(1) + r(2) * r(2) + r(3) * r(3));
	rdot(1) = r(4);
	rdot(2) = r(5);
	rdot(3) = r(6);
	rdot(4) = - (GRAVITY_CONST / delta) * (EARTH_MASS / delta) * (r(1) / delta);
	rdot(5) = - (GRAVITY_CONST / delta) * (EARTH_MASS / delta) * (r(2) / delta);
	rdot(6) = - (GRAVITY_CONST / delta) * (EARTH_MASS / delta) * (r(3) / delta);
endfunction
