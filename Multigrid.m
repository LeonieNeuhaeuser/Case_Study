% Assume A to be square and  its dimension is odd 
function [u1_m_grid, errvect,iter]=Multigrid(A,F,u0,tol)

    % Initialization
    err = norm(F-A*u0,2);
    errvect(1) = err;
    k = 1;
    u_i = u0;
    [R, P, A_bar] = multigrid_operator(A);
    n = size(A);
    steps = 5;

    while err >  tol
        % Pre-smoothing
        u_s = smoothing(A,F,u_i,steps,1/2);

        % Calculate  residual in the fine grid 
        r_s = F - A*u_s;

        % Restrict residual to coarse grid 
        r_s_bar = R*r_s;

        %Calculatre the correction
        e_s_bar = A_bar\r_s_bar;

        % Prolong residual to fine grid
        e_s = P*e_s_bar;

        %Update solution
        u_i = u_s + e_s;

        % Post-smoothing 
        u_i = smoothing(A,F,u_i,steps,1/2);

        err = norm(F-A*u_i,2);
        k = k+1;
        errvect(k) = err;

    end

    iter =k-1;
    u1_m_grid =u_i;
end 