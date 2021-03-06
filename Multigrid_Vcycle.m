% Assume A to be square and  its dimension is odd 
function [u, errvect] = Multigrid_Vcycle(A,F,u0,tol,presteps,poststeps,nmin,level,max_depth)

    u = u0;
    err0 = norm(A*u-F,2);
    err = err0;
    errvect = [1];
    [R, P, A_bar] = multigrid_operator(A);
    coarse_n = size(R,1);
    n = size(F,1);

    if n <= nmin || level == max_depth
        u = A\F;
        return;
    end
    
    while err > tol
            
    % Pre-smoothing
    u_s = GS_smoothing(A,F,u,presteps,1);

    % Calculate  residual in the fine grid 
    r_s = F - A*u_s;

    % Restrict residual to coarse grid 
    r_s_bar = R*r_s;

    %Solve coarse grid with recursion
    v0 = zeros(coarse_n,1);
    e_s_bar = Multigrid_Vcycle(A_bar, r_s_bar, v0, tol, presteps, poststeps,nmin, level + 1, max_depth);

    % Prolong residual to fine grid
    e_s = P*e_s_bar;

    %Update solution
    u = u_s + e_s;
    
    % Post-smoothing 
    u = GS_smoothing(A,F,u,poststeps,1);
    
    %error 
    err = norm(A*u-F,2)/err0;
    errvect = [errvect err];
    end
    
    errvect = err0*errvect; 

end
        
