%~ Returns the Rot matrix that transforms coordinates from the orbital plane to the 
%~ reference plane.


function Rot = calc_rotmatrix (sinw, cosw, i, omega)

	P_x = cosw * cos(omega) - sinw * sin(omega) * cos(i);
	P_y = cosw * sin(omega) + sinw * cos(omega) * cos(i);
	P_z = sinw * sin(i);
	Q_x = -sinw * cos(omega) - cosw * sin(omega) * cos(i);
	Q_y = -sinw * sin(omega) + cosw * cos(omega) * cos(i);
	Q_z = cosw * sin(i);

	Rot = [...
			P_x, Q_x; ...
			P_y, Q_y; ...
			P_z, Q_z ...
			]';
end

