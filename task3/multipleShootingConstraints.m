function [c, ceq] = multipleShootingConstraints(z, U0_idx, Uj_idx, Y_inds, N, nx, Y0_INIT, dt)
    ceq = [];

    % Initial node equality: Y0 == Y0_INIT
    Y0 = z( Y_inds(0) );
    ceq = [ceq; Y0 - Y0_INIT];

    % Continuity constraints for each interval: Y_{k+1} - Phi(Y_k, u_k) = 0
    for k = 1:N
        Yk  = z( Y_inds(k-1) );
        Yk1 = z( Y_inds(k) );
        uk  = [ z(U0_idx(k)); z(Uj_idx(k)) ];
        PhiY = solveEquationSystemTask3(Yk, uk, dt);   % uses daeSystemLHS(t,y,PARAMS)
        ceq = [ceq; Yk1 - PhiY];
    end

    % Rate limits on controls: |ΔU| <= 5 K between consecutive intervals
    c = [];
    if N >= 2
        dU0 = zeros(N-1,1);
        dUj = zeros(N-1,1);
        for k = 1:N-1
            dU0(k) = z(U0_idx(k+1)) - z(U0_idx(k));
            dUj(k) = z(Uj_idx(k+1)) - z(Uj_idx(k));
        end
        % Convert two-sided to <= form
        c = [ c;
              dU0 - 5;         %  dU0 <= 5
             -dU0 - 5;         % -dU0 <= 5  (i.e., dU0 >= -5)
              dUj - 5;
             -dUj - 5 ];
    end
end