% Satellite trajectory plot

function [] = plot_orbit(xyz, range)
    xyz = xyz / 1e6;
    plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3));
    axis(range);
    view(3);
    xlabel('x [1000 km]');
    ylabel('y [1000 km]');
    zlabel('z [1000 km]');
    legend('location', 'northeastoutside');
end
