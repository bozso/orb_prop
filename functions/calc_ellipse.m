% Returns the theoretical coordinates (x,y) in the orbital plane.

% Inputs: orbital elements: tau, e, a, n
% Outputs: x, y coordinate pairs in the orbital plane in matrix form (ellipse)

function ellipse = calc_ellipse (time, tau, e, a, n)
    M = n * (time - tau);

    for iii = 1 : numel(time)
            [E(iii), sinv(iii), cosv(iii)] = newtonm_scv (e, M(iii));
    end
    
    r = a * (1 - e * cos(E));
    ellipse(:,1) = r .* cosv;
    ellipse(:,2) = r .* sinv;
end
