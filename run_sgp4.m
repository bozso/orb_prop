	t_r_ode = load('-ascii', 'output/4793.dat');
	
	time = t_r_ode(:,1);
	
	clear t_r_ode;
	
	outfile = fopen('output/4793_sgp4.dat', 'w');
	
	infile = fopen('/home/istvan/orb_prop/synthetic.tle', 'r');
	
	fgetl(infile);
	
	longstr1 = fgets(infile, 130);
	while ( (longstr1(1) == '#') && (feof(infile) == 0) )
		longstr1 = fgets(infile, 130);
	end
	
	if (feof(infile) == 0)
		longstr2 = fgets(infile, 130);
	end
		
	satrec = process(whichconst, longstr1, longstr2);
	
	fclose(infile);
	
	global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  
	
	global opsmode

	for iii = 1:numel(time)
		[satrec, ro ,vo] = sgp4 (satrec,  time(iii));

		if (satrec.error > 0)
			fprintf(1,'# *** error: t:= %f *** code = %3i\n', time(iii), satrec.error);
		end  
	                
		if (satrec.error == 0)
			fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n',...
					satrec.t,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));
		else
			jd = satrec.jdsatepoch + time(iii)/1440.0;
			[year,mon,day,hr,minute,sec] = invjday ( jd );
			
			fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f',...
					tsince,ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));
			
			[p,a,ecc,incl,node,argp,nu,m,arglat,truelon,lonper ] = rv2coe (ro,vo,mu);
	
			fprintf(outfile, ' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f %5i%3i%3i %2i:%2i:%9.6f \n',...
				a, ecc, incl*rad, node*rad, argp*rad, nu*rad, m*rad,year,mon,day,hr,minute,sec );
		end %// if satrec.error == 0

	end %// while propagating the orbit   
	
	fclose(outfile);
