function [t_all, Y_all] = simulate_trajectory(Y0, U0_seq, Uj_seq, dt)
    N = numel(U0_seq);
    MASS = eye(4);
    opts = odeset('Mass', MASS, 'RelTol', 1e-7, 'AbsTol', 1e-8);

    t_all = [];
    Y_all = [];
    Yk = Y0;
    t_offset = 0;

    for k = 1:N
        PARAMS = [U0_seq(k), Uj_seq(k)];
        odefun = @(t, y) daeSystemLHS(t, y, PARAMS);
        [t, Y] = ode15s(odefun, [0 dt], Yk, opts);

        if k == 1
            t_all = [t_all; t];
            Y_all = [Y_all; Y];
        else
            t_all = [t_all; t_offset + t(2:end)];
            Y_all = [Y_all; Y(2:end,:)];
        end

        t_offset = t_offset + dt;
        Yk = Y(end,:).';
    end
end