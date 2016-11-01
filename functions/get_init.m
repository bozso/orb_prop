function rv0 = get_init(satrec)

    % Getting initial conditions in TEME
    printf("Getting initial conditions from satrec...");
    [satrec, rteme ,vteme] = sgp4(satrec,  0.0);

    % Transforming into ECI
    t_poz_vel_0 = teme2eci_full([satrec.jdsatepoch, rteme ,vteme]);

    % Changing into SI units
    rv0 = [t_poz_vel_0(2:4) * 1e3, t_poz_vel_0(5:end) * 1e3];

end




