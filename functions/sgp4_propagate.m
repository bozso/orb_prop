% Satellite trajectory calculation with SGP

function [t_poz_vel] = sgp4_propagate(step, day, model, satrec)
    global whichconst day2sec sec2day
	
	time = 0:step:day*day2sec;
	
	t_poz_vel = zeros(numel(time), 7);
	
    switch (model)
        case 'point_mass'
            printf("sgp4_propagate: SGP4 can not handle the Earth as ");
            printf("a point mass.\n");
%            return
        case 'zonal'
            printf("sgp4_propagate: Zonal harmonics will be used, air ");
            printf("friction will be inored.\n");
            satrec.bstar = 0.0;
        otherwise
            disp('sgp4_propagate: Unrecognized model option!');
            return
    end
    
    global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  
    
    global opsmode

	sec2min = 1 / 60;

    for iii = 1:numel(time)
        [satrec, ro ,vo] = sgp4 (satrec,  time(iii) * sec2min);

        if (satrec.error > 0)
            fprintf(1,'# *** error: t:= %f *** code = %3i\n', ...
            time(iii), satrec.error);
        end  
                
        if (satrec.error == 0)
            t_poz_vel(iii,:) = [satrec.jdsatepoch + time(iii) * sec2day,...
                                ro, vo];
%            [satrec.jdsatepoch + time(iii) * sec2day, ro, vo]
%            t_poz_vel(iii,:)
        else
            printf("sgp4_propagate: Error: satrec.error is non zero! ");
            printf("Exiting now!\n");
            return 
        end %// if satrec.error == 0
    end %// while propagating the orbit   
end	
