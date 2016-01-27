sgp4 = load('-ascii', 'output/sgp4_fengyun.dat');
rk = load('-ascii', 'output/25730.dat');

delta = [rk(:,2:4) / 1e3 - sgp4(:,2:4)];

delta = sqrt(sum(delta.^2, 2));

h = figure(1);
xlabel('t [day]')
ylabel('delta [m]')
plot(rk(:,1) / 86400, delta)
print(h, '-dpng', 'sgp4_vs_rk.png');
