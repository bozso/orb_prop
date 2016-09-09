function [jday_TT] = jday2jday_tt(jday_UT1)
	global eod deltaT_jday
		
	TT_minus_UTC = interp1(deltaT_jday(:,1), deltaT_jday(:,2), ...
	jday_UT1) + interp1(eod(:,1), eod(:,6), jday_UT1);
	
	jday_TT = jday_UT1 + TT_minus_UTC ./ 86400;
end
