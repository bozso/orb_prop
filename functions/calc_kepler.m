function [r, v] = calc_kepler(r0, v0, step, nsteps)

    r = zeros(nsteps, 3);
    v = zeros(nsteps, 3);

    r(1,:) = r0;
    v(1,:) = v0;

    for iii = 2:nsteps
        [r(iii,:), v(iii,:)] =  kepler(r(iii - 1,:) , v(iii - 1,:), step);
    end

end
