% Satellite trajectory plot

function [] = plot_orbit(xyz, figure_num)
    h = figure(figure_num)
    plot3(xyz(:, 1) / 1e6, xyz(:, 2) / 1e6, xyz(:, 3) / 1e6, 'ob', ...
            'linewidth', 1, 0, 0, 0, 'ok;Earth;');
    view(3);
    xlabel('x [1000 km]');
    ylabel('y [1000 km]');
    zlabel('z [1000 km]');
    legend('location', 'northeastoutside');
end
