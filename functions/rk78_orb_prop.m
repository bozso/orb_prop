%~ Integrates the equations of motion for given time intervals

function timexyzv = rk78_orb_prop (r_0, timestep, cicles, vopt)
	
	[t, r_rk4] = ode78(@Gravi, [0, timestep], r_0, vopt);
	
	xyzv = r_rk4;
	time = t;
	
	r_0 = r_rk4(rows(r_rk4), :);
	
	for iii = 2:cicles
		[t, r_rk4] = ode78(@Gravi, [(iii - 1) * timestep, iii * timestep], r_0, vopt);
	
		xyzv = [xyzv; r_rk4(2:end, :)];
		time = [time; t(2:end, :)];

		r_0 = r_rk4(end, :);
	end

	timexyzv = [time, xyzv];
end

