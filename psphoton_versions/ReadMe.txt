Version 0 :
	v0 is the original version from Mancinelli, using Earth radius. Directly form GitHub.
Version 1 (Update Based on V0) :
	v1 changes the radius 'erad' to 3390;
	v1 changes the PSDF from Von Karman type to Henry-Green type. Modified in function 'expsato'
	v1 adds some print commands to better display the process.
	v1 can output the flattened velocities and depth to output file.
Version 1.1 :
	v1.1 increases the integration points in function GAVGSATO2 from n=60 to n=400 to have 		finer azimuth and elevation calculation
Version 1.2 :
	v1.2 Takes maximum rms into ray parameter calculation : pmax = 1./(vsmin * (1-rms_max)) 
Version 1.3 (Update based on V1.2)
	v1.3 Changes the intrinsic quality factor calculation to Qp = 5Qs 
Version 1.3.1 (Update based on V1.3)
	v1.3.1 makes the ratio Qs/Qp as an input parameter. And add safety check
Version 1.3.2 (Update based on V1.3.1)
	v1.3.2 removes the single input Qs2Qp, but adds new input parameter: Qs2Qp for each attenuation layer.
Version 1.3.3 (Update based on V1.3.2)
	v1.3.3 adds Hurst exponent as input for each scattering layer.
Version 1.3.4 (Update based on V1.3.3)
	v1.3.4 fixes the PSDF function to become gamma(k+1.5), not gamma(k+0.5)


Version 2 (Update based on V1) :
	v2 adds options to write ray parameter output 'iwrite_rayparam'
	v2 adds option to write scattering event output 'iwrite_scatevent
Version 2.1 (Update based on V2): 
	v2.1 increases the integration points in function GAVGSATO2 from n=60 to n=400 to have 		finer azimuth and elevation calculation
Version 2.2 (Update based on V2.1) :
	v2.2 changes the slowness calculation by taking the maximum rms into calculation
	v2.2 adds print information when the Vs=0 at scatterers
	v2.2 adds error log when 'SCATRAYPOL problema' happens
