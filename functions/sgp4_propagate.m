function [] = sgp4_propagate(infile_path, outname)
    global whichconst
    
	[step day model satrec] = process_infile(infile_path, whichconst);
	
	time = 0:step:day*86400;
	
    switch (model)
            case 'point_mass'
                    disp('sgp4_propagate: SGP4 can not handle the Earth as a point mass.');
                    return
            case 'zonal'
                    disp('sgp4_propagate: Zonal harmonics will be used, air friction will be inored')
                    satrec.bstar = 0.0;
            otherwise
                    disp('sgp4_propagate: Unrecognized model option!');
                    return
    end
    
    global tumin mu radiusearthkm xke j2 j3 j4 j3oj2  
    
    global opsmode

    outfile = fopen(['output/', outname, '_sgp.dat'], 'w');

    for iii = 1:numel(time)
            [satrec, ro ,vo] = sgp4 (satrec,  time(iii) / 60.0);

            if (satrec.error > 0)
                    fprintf(1,'# *** error: t:= %f *** code = %3i\n', time(iii), satrec.error);
            end  
                    
            if (satrec.error == 0)
                    fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n',...
                                    time(iii),ro(1)*1e3,ro(2)*1e3,ro(3)*1e3,vo(1),vo(2),vo(3));
            else
                    jd = satrec.jdsatepoch + time(iii)/1440.0;
                    [year,mon,day,hr,minute,sec] = invjday ( jd );
                    
                    fprintf(outfile, ' %16.8f %16.8f %16.8f %16.8f %12.9f %12.9f %12.9f',...
                                    time(iii),ro(1)*1e3,ro(2)*1e3,ro(3)*1e3,vo(1),vo(2),vo(3));
                    
                    [p,a,ecc,incl,node,argp,nu,m,arglat,truelon,lonper ] = rv2coe (ro,vo,mu);
    
                    fprintf(outfile, ' %14.6f %8.6f %10.5f %10.5f %10.5f %10.5f %10.5f %5i%3i%3i %2i:%2i:%9.6f \n',...
                            a, ecc, incl*rad, node*rad, argp*rad, nu*rad, m*rad,year,mon,day,hr,minute,sec );
            end %// if satrec.error == 0

    end %// while propagating the orbit   
    
    fclose(outfile);
end	
	
	
