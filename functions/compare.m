% Compares two satellite trajectories

function [] = compare(t_r_1, t_r_2, outname, figure_num)
    if ( rows(t_r_1) != rows(t_r_2) )
        disp('compare: Rows of the two data matrices are not equal!')
        return
    end

    if ( sum( t_r_1(:,1) - t_r_2(:,1) ) ~= 0 )
        disp('compare: Dates of the two data matrices are not equal!')
        return
    end

    delta = [t_r_1(:,2:4) - t_r_2(:,2:4)];
    delta = sqrt(sum(delta.^2, 2));
    delta = [t_r_1(:,1) delta];
    save('-ascii', ['output/', outname], 'delta');

	figure(figure_num)
    hold on
    plot(delta(:,1) / 86400, delta(:,2), '-')
    xlabel('t [day]')
    ylabel('delta [m]')
    hold off
end

