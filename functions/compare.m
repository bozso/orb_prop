function [] = compare(t_r_1, t_r_2, outname)
    if ( rows(t_r_1) != rows(t_r_2) )
        disp('Rows of the two data matrices are not equal!')
        return 0
    end

    if ( sum( t_r_1(:,1) == t_r_2(:,2) ) != 0 )
        disp('Dates of the two data matrices are not equal!')
        return 0
    end

    delta = [t_r_1(:,2:4) - t_r_2(:,2:4)];
    delta = sqrt(sum(delta.^2, 2));
    delta = [t_r_1(:,1) delta];
    save('-ascii', outname, 'delta');

    hold on
    plot(delta(:,1) / 86400, delta(:,2), '-')
    xlabel('t [day]')
    ylabel('delta [m]')
    hold off
end

