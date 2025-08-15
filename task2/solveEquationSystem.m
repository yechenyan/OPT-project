function[X,Y]=solveEquationSystem(PARAMS)
    %% Initial conditions at t=0
    CA_init = 10;    % mol/m^3
    CB_init = 18;    % mol/m^3
    CC_init = 0;     % mol/m^3
    T_init = 343;     % K
    Y_INIT = [CA_init; CB_init; CC_init; T_init];

    
    %% Time span
    X_INTERVAL = [0 3600];  % seconds

    MASS = eye(4,4);
    OPTIONS = odeset('Mass',MASS, 'RelTol',1e-7,'AbsTol',1e-9);

    odefun = @(t, y) daeSystemLHS(t, y, PARAMS);
    [X, Y] = ode15s(odefun, X_INTERVAL, Y_INIT, OPTIONS);

end