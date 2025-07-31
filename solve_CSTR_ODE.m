

function dydt = solve_CSTR_ODE(~, y, T0, Tj)
   
    %% Constants
    EA = 69.1e3;     % J/mol
    EB = 75.2e3;     % J/mol
    kA0 = 5.0e6;     % 1/s
    kB0 = 1.6e6;     % 1/s
    rho = 800;       % kg/m^3
    cp = 3.5e3;      % J/kg/K
    UA = 1.4;        % W/K
    deltaHA = -45.1e3; % J/mol
    deltaHB = -55.1e3; % J/mol
    F = 3e-3;        % m^3/s
    V = 3e-3;        % m^3

    %% Feed conditions
    CA0 = 5;         % mol/m^3
    CB0 = 15;        % mol/m^3

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
    dcBdt = (F/V)*(CB0 - cB) - kB * cB + kA * cA;
    dcCdt = -(F/V) * cC + kA * cA;

    % Energy balance
    dTdt = (F/V)*(T0 - T) ...
           + (UA / (rho * cp * V)) * (Tj - T) ...
           + (-deltaHA * kA * cA + (-deltaHB) * kB * cB) / (rho * cp);

    % Pack the derivatives
    dydt = [dcAdt; dcBdt; dcCdt; dTdt];
end