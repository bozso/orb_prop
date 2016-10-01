% Processes inputfile

function [step day model satrec] = process_infile(file_loc)
    global whichconst jday_low jday_high
    
    % File that contains propagation parameters and initial conditions
    infile = fopen(file_loc, 'r');
    
    % Processing parameters and initial conditions
    
    str = strsplit(fgetl(infile));
    
    if (numel(str) > 3)
            disp('Error: Too many options given in the first line!')
            return
    end
    
    step = str2num(str{1});
    day = str2num(str{2});
    model = str{3};
    
    % Processing oribtal element set
    printf('Reading satellite data...')
    longstr1 = fgets(infile, 130);
    while ( (longstr1(1) == '#') && (feof(infile) == 0) )
            longstr1 = fgets(infile, 130);
    end
    
    if (feof(infile) == 0)
            longstr2 = fgets(infile, 130);
    end
            
#    satrec = process(whichconst, longstr1, longstr2);
#    [satrec, startmfe, stopmfe, deltamin] = twoline2rv(whichconst, longstr1, ...
#          longstr2, 'c','e');
    [startmfe, stopmfe, deltamin, satrec] = twoline2rv(longstr1, longstr2, ...
          'c', 'm', 'i', whichconst);
    fclose(infile);
    printf('DONE\n');

    if (satrec.jdsatepoch > jday_high || satrec.jdsatepoch < jday_low)
	{
        printf("process_infile: Warning: Can only compare sgp and "); 
        printf("Runge-Kutta if ephemeris date is between 1992.01.01 ");
        printf("2016.07.01\n");
	}	
end
