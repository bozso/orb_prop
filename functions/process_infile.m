function [step day model satrec] = process_infile(file_loc, whichconst)
	%~ File that contains propagation parameters and initial conditions
	infile = fopen(file_loc, 'r');
	
	%~ Processing parameters and initial conditions
	
	str = strsplit(fgetl(infile));
	
	if (numel(str) > 3)
		disp('Error: Too many options given in the first line!')
		return
	end
	
	step = str2num(str{1});
	day = str2num(str{2});
	model = str{3};
	
	%~ Processing oribtal element set
	disp('Reading satellite data...')
	longstr1 = fgets(infile, 130);
	while ( (longstr1(1) == '#') && (feof(infile) == 0) )
		longstr1 = fgets(infile, 130);
	end
	
	if (feof(infile) == 0)
		longstr2 = fgets(infile, 130);
	end
		
	satrec = process(whichconst, longstr1, longstr2);
	
	fclose(infile);
end

