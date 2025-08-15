
function Yend = solveEquationSystemTask3(Ystart, u, dt)
    PARAMS = [u(1), u(2)];
    odefun = @(t, y) daeSystemLHS(t, y, PARAMS);  % <-- use original model signature
    MASS = eye(4);
    opts = odeset('Mass', MASS, 'RelTol', 1e-7, 'AbsTol', 1e-8);
    [~, Y] = ode15s(odefun, [0 dt], Ystart, opts);
    Yend = Y(end,:).';
end