% Returns the Rot matrix that transforms coordinates from the orbital plane to the 
% reference plane.


function Rot = calc_rotmatrix (sinw, cosw, incl, omega)
    P_x = cosw * cos(omega) - sinw * sin(omega) * cos(incl);
    P_y = cosw * sin(omega) + sinw * cos(omega) * cos(incl);
    P_z = sinw * sin(incl);
    Q_x = -sinw * cos(omega) - cosw * sin(omega) * cos(incl);
    Q_y = -sinw * sin(omega) + cosw * cos(omega) * cos(incl);
    Q_z = cosw * sin(incl);

    Rot =  [P_x, Q_x; ...
            P_y, Q_y; ...
            P_z, Q_z]';
end

