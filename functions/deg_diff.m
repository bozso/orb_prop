function delta_sigma = deg_diff(r1, r2)

    row = rows(r1);

    delta_sigma = zeros(row, 1);

    [theta_1, phi_1, r_1] = cart2sph(r1);
    [theta_2, phi_2, r_2] = cart2sph(r2);

    for iii = 1:rows(r1)
        delta_sigma(iii) = polar_dist(theta_1(iii), theta_2(iii), ...
        phi_1(iii), phi_2(iii));
    end

end
