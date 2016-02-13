% Equations of motion, assuming that the Earth is a point mass

function rdot = point_mass(t, r)
    global mu_si
    delta = norm([r(1), r(2), r(3)]);
    mu_per_rcubed = mu_si / delta^3;
    rdot(1) = r(4);
    rdot(2) = r(5);
    rdot(3) = r(6);
    rdot(4) = - mu_per_rcubed * r(1);
    rdot(5) = - mu_per_rcubed * r(2);
    rdot(6) = - mu_per_rcubed * r(3);
end
