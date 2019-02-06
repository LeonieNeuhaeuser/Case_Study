function [u ,err,errvec]= GS(A,f,u0,tol)

    n = length(u0);
    u = u0; 
    err = max(abs(A*u-f));
    temp = zeros(n,1);
    diagonal = diag(A);
    errvec =[];
    
    while err > tol
        for i = 1:n                        
            u(i) = 1/diagonal(i)*(-dot(A(i,1:i-1),u(1:i-1))-dot(A(i,i+1:n),u(i+1:n))+f(i));
        end  
        err = max(abs(A*u-f));        
        errvec = [errvec err];
    end
    

end