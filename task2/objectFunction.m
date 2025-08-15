function c_a_end = objectFunction (params)
    [~, Y] = solveEquationSystem(params);
    c_a_end = -Y(end, 1);   % maximize CA(tf)
end