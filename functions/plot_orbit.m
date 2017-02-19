% Satellite trajectory plot

function [] = plot_orbit(xyz, range, unit)
    plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3));

    if (numel(range) > 1)
        axis(range);
    end
    view(3);
    xlabel(['x [' unit ']']);
    ylabel(['y [' unit ']']);
    zlabel(['z [' unit ']']);
end
