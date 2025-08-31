function dydt = daeSystemLHS(X, y, PARAMS)

    T0 = PARAMS(1);
    Tj = PARAMS(2);
   
    %% Constants
    EA      = 69000.0;    % J/mol
    EB      = 72000.0;    % J/mol
    kA0     = 5.0e6;      % 1/s
    kB0     = 1.0e7;      % 1/s
    rho     = 800.0;      % kg/m^3
    cp      = 3.5;        % J/kg/K 
    UA      = 1.4;        % W/K
    deltaHA =  45000.0;   % J/mol   
    deltaHB = -55000.0;   % J/mol
    F       = 6.5e-4;     % m^3/s
    V       = 1.0;        % m^3
    
    %% Feed conditions
    CA0 = 5.0;            % mol/m^3
    CB0 = 15.0;           % mol/m^3

    % Unpack state variables
    cA = y(1);
    cB = y(2);
    cC = y(3);
    T  = y(4);

    % Rate constants with Arrhenius law
    kA = kA0 * exp(-EA / (8.314 * T));
    kB = kB0 * exp(-EB / (8.314 * T));

    % Material balances
    dcAdt = (F/V)*(CA0 - cA) - kA * cA + kB * cB;
    dcBdt = (F/V)*(CB0 - cB) - kB * cB; 
    dcCdt = -(F/V) * cC + kA * cA;

    % Energy balance
    dTdt = (F/V)*(T0 - T) ...
           + (UA / (rho * cp * V)) * (Tj - T) ...
           + (-deltaHA * kA * cA + (-deltaHB) * kB * cB) / (rho * cp);

    % Pack the derivatives
    dydt = [dcAdt; dcBdt; dcCdt; dTdt];
end