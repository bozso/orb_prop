function [t_poz_vel_eci] = teme2eci(t_poz_vel_teme)
	global eod deltaT_jday
	
	jday_UT1 = t_poz_vel_teme(:,1);
	
	ttt = jday2jday_tt(jday_UT1) ./ 36525;
	
	lod = interp1(eod(:,1), eod(:,8), jday_UT1) .* 1e-3;
	xp = interp1(eod(:,1), eod(:,2), jday_UT1) ./ 206265;
	yp = interp1(eod(:,1), eod(:,4), jday_UT1) ./ 206265;
	
	num_of_rows = size(t_poz_vel_teme, 1);
	
	recef = zeros(3, num_of_rows);
	vecef = zeros(3, num_of_rows);
	
	for iii = 1:num_of_rows
		[r_ecef(:,iii), v_ecef(:,iii)] = teme2ecef(...
											t_poz_vel_teme(iii,2:4)',...
											t_poz_vel_teme(iii,5:7)',...
											ttt(iii), jday_UT1(iii), lod(iii),...
											xp(iii), yp(iii));
	end
	
	
	[r_eci v_eci] = ECEFtoECI(jday_UT1, r_ecef, v_ecef);
	
	t_poz_vel_eci = [jday_UT1, r_eci', v_eci'];
end
