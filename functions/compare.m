% Compares two satellite trajectories

function [delta] = compare(poz_1, poz_2)
    if ( size(poz_1) ~= size(poz_2) )
        disp('compare: Size of the two data matrices are not equal!')
        return
    end
    
    delta = poz_1 - poz_2;
    delta = sqrt(sum(delta.^2, 2));
end

