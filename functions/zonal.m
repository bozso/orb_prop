function rdot = zonal(t, r)
    r_pos = [r(1) r(2) r(3)];
    rdot(1) = r(4);
    rdot(2) = r(5);
    rdot(3) = r(6);
    
    h = 1;
    
    rdot(4) = (pot_zonal(r_pos + [h, 0, 0]) - pot_zonal(r_pos - [h, 0, 0])) / (2 * h);
    rdot(5) = (pot_zonal(r_pos + [0, h, 0]) - pot_zonal(r_pos - [0, h, 0])) / (2 * h);
    rdot(6) = (pot_zonal(r_pos + [0, 0, h]) - pot_zonal(r_pos - [0, 0, h])) / (2 * h);
end
