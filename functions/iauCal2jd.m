%--------------------------------------------------------------------------
% Purpose:
%
%   Modified Julian Date from calendar date and time
%
% Input(s):
%
%   Year      
%   Month
%   Day
%   Hour      
%   Min
%   Sec
% 
% output(s):  
%    
%    djm0  
%    djm     Modified Julian Date
% 
% Reference: O. Montenbruck, E. Gill (2000), "Satellite Orbits"
% 
%--------------------------------------------------------------------------

function [djm0 djm] = iauCal2jd(Year, Month, Day, Hour, Min, Sec)

  if (nargin == 3)
      Hour = 0;
      Min = 0;
      Sec = 0;
  end

  djm0 = 2400000.5;
   
  if ( Month<=2 ) 
      Month = Month+12;
      Year = Year-1;
  end
  
  if ( (10000*Year+100*Month+Day) <= 15821004 )
    b = fix(-2 + ((Year+4716)/4) - 1179);     % Julian calendar
  else
    b = fix((Year/400)-(Year/100)+(Year/4));  % Gregorian calendar
  end
  
  MjdMidnight = 365*Year - 679004 + b + fix(30.6001*(Month+1)) + Day;
  FracOfDay   = (Hour+Min/60.0+Sec/3600.0) / 24.0;
  
  djm = MjdMidnight + FracOfDay;
