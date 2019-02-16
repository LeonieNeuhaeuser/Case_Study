% Assume A to be square and  its dimension is odd 
function [u, errvect, k] = Generalised_multigrid(A,F,u0,tol,presteps,poststeps,level,levels)

    % Initialisation
    err = norm(F-A*u0,2);
    errvect = [err];
    k = 0;
    u = u0;
    
    while err > tol 
       
        % Run one Vcycle
        u = Multigrid_Vcycle(A,F,u,presteps,poststeps,level,levels);
        
        % Calculate error
        err = norm(F-A*u,2); 
        errvect = [errvect err]; 
        k = k + 1;   
        
end 