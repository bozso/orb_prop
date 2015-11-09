function rdot = zonal(t, r)
	rdot(1) = r(4);
	rdot(2) = r(5);
	rdot(3) = r(6);
	
	h = 5;
	
	g(1) = (pot_zonal(r + [h, 0, 0]) - pot_zonal(r - [h, 0, 0])) / (2 * h);
	g(2) = (pot_zonal(r + [0, h, 0]) - pot_zonal(r - [0, h, 0])) / (2 * h);
	g(3) = (pot_zonal(r + [0, 0, h]) - pot_zonal(r - [0, 0, h])) / (2 * h);
	
	rdot(4) = g(1);
	rdot(5) = g(2);
	rdot(6) = g(3);
end
