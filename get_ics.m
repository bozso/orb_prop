addpath ('/home/istvan/orb_prop/sgp4');

global tumin radiusearthkm xke j2 j3 j4 j3oj2
[tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2] = getgravc(84);

infile = fopen(argv(){1});

longstr1 = fgets(infile, 130);
while ( (longstr1(1) == '#') && (feof(infile) == 0) )
	longstr1 = fgets(infile, 130);
end

if (feof(infile) == 0)
	longstr2 = fgets(infile, 130);
end

satrec = process(84, longstr1, longstr2);

%~ call the propagator to get the initial state vector value
[satrec, ro ,vo] = sgp4 (satrec,  0.0);

ro = ro * 1e3;
vo = vo * 1e3;

printf('%16.8f %16.8f %16.8f %12.9f %12.9f %12.9f\n',...
                 ro(1),ro(2),ro(3),vo(1),vo(2),vo(3));
