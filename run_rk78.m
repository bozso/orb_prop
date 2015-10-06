GRAVITY_CONST = 6.67384e-11;
EARTH_MASS = 5.97237e24;

addpath('/home/istvan/orb_prop_matsource/');
addpath('/home/istvan/orb_prop_matsource/rk_orbit/functions');

%r_0 = [4e6; 5.74456e6; 0; -4571.43; 6565.214; 0]';
%r_0 = [3.5e6; ; 3.5e6; 4949770; 0; 8e3; 0]';
r_0 = [7e6; 0; 0; 0; 8e3; 0]';

plotname = sprintf('output/v_%.2f_%.2f_%.2f', r_0(4), r_0(5), r_0(6));

vopt = odeset ('RelTol', 1e-4, 'AbsTol', 1, 'NormControl', 'on', 'InitialStep', 125, 'MaxStep', 200);

[a, e, inclination, omega, sinw, cosw, n, tau] = OrbitalElementSet(r_0, 0.0, GRAVITY_CONST * EARTH_MASS);

Rot = CalcRotMtx (sinw, cosw, inclination, omega);

step = 100;

day = 4;

cicles = ceil((86400 * day ) / step);

timexyzv = rk78_orb_prop (r_0, step, cicles, vopt);

xyz = [timexyzv(:,2), timexyzv(:,3), timexyzv(:,4)];

time = timexyzv(:, 1);

save('-ascii', [plotname, '.dat'], 'timexyzv');

clear timexyzv;

ellipse = CalcEllipse (time, tau, e, a, n);

ellipse = ellipse * Rot;

save ('-ascii', [plotname, '_ellipse.dat'], 'ellipse');

delta = ellipse - xyz;

delta = sqrt(sum(delta'.^2));

delta = [time, delta'];

h=figure(1);
hold on;
plot3 (xyz(:, 1) / 1e6, xyz(:, 2) / 1e6, xyz(:, 3) / 1e6, 'ob;RK4;', 'linewidth', 1, ...
ellipse(:,1) / 1e6, ellipse(:,2) / 1e6, ellipse(:,3) / 1e6, '*r;Ellipse;', 'linewidth', 1, ...
0, 0, 0, 'ok;Earth;');
view(3);
xlabel('x [1000 km]');
ylabel('y [1000 km]');
zlabel('z [1000 km]');
legend ('location', 'northeastoutside');
print(h,'-dpng','-color', [plotname, '.png']);
hold off;

save('-ascii', [plotname, '_delta.dat'], 'delta');

h = figure(2);
clf;
hold on
plot (delta(:,1) / 86400, delta(:,2));
print(h,'-dpng','-color', [plotname, '_delta' , '.png']);
hold off

