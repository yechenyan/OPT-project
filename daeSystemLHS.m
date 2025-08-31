function DYDX = daeSystemLHS(X, Y, VARS)
   
    % declare parameters 
    PARAMS(1) = 15.0;  	% e0_CB0 
	PARAMS(2) = 72000.0;  	% e0_EB 
	PARAMS(3) = 6.5E-4;  	% e0_F 
	PARAMS(4) = 1.0;  	% e0_V 
	PARAMS(5) = 1.0E7;  	% e0_kB0 
	PARAMS(6) = 69000.0;  	% e0_EA 
	PARAMS(7) = 5000000.0;  	% e0_kA0 
	PARAMS(8) = 5.0;  	% e0_CA0 
	PARAMS(9) = VARS(1);  	% e0_T0 
	PARAMS(10) = VARS(2);  	% e0_Tj 
	PARAMS(11) = 1.4;  	% e0_UA 
	PARAMS(12) = 3.5;  	% e0_cp 
	PARAMS(13) = 45000.0;  	% e0_deltaHA 
	PARAMS(14) = -55000.0;  	% e0_deltaHB 
	PARAMS(15) = 800.0;  	% e0_rho 

	% read out variables  
	 
    e0_cA = Y(1); 
	e0_cB = Y(2); 
	e0_cC = Y(3); 
    e0_T = Y(4);

	% read out differential variable
	e0_t = X;

	% read out parameters  
	e0_CB0 = PARAMS(1); 
	e0_EB = PARAMS(2); 
	e0_F = PARAMS(3); 
	e0_V = PARAMS(4); 
	e0_kB0 = PARAMS(5); 
	e0_EA = PARAMS(6); 
	e0_kA0 = PARAMS(7); 
	e0_CA0 = PARAMS(8); 
	e0_T0 = PARAMS(9); 
	e0_Tj = PARAMS(10); 
	e0_UA = PARAMS(11); 
	e0_cp = PARAMS(12); 
	e0_deltaHA = PARAMS(13); 
	e0_deltaHB = PARAMS(14); 
	e0_rho = PARAMS(15);

	% evaluate the function values  
	DYDX(1) = ((e0_F)/(e0_V)) * (e0_T0 - e0_T) + ((e0_UA)/((e0_rho * e0_cp * e0_V))) * (e0_Tj - e0_T) + (( - e0_deltaHA * e0_kA0 * exp( - (e0_EA)/((8.314 * e0_T))) * e0_cA + ( - e0_deltaHB) * e0_kB0 * exp( - (e0_EB)/((8.314 * e0_T))) * e0_cB))/((e0_rho * e0_cp)); 
	DYDX(2) = ((e0_F)/(e0_V)) * (e0_CB0 - e0_cB) - e0_kB0 * exp( - (e0_EB)/((8.314 * e0_T))) * e0_cB; 
	DYDX(3) = ((e0_F)/(e0_V)) * (e0_CA0 - e0_cA) - e0_kA0 * exp( - (e0_EA)/((8.314 * e0_T))) * e0_cA + e0_kB0 * exp( - (e0_EB)/((8.314 * e0_T))) * e0_cB; 
	DYDX(4) =  - ((e0_F)/(e0_V)) * e0_cC + e0_kA0 * exp( - (e0_EA)/((8.314 * e0_T))) * e0_cA; 

	DYDX=DYDX';
end