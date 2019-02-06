function [u ,err]= SSOR(w,A,f,u0,tol)

    n = length(u0);
    u = u0; 
    err = max(abs(A*u-f));
    temp = zeros(n,1);
    
    while err > tol
        for i = 1:n            
            a = A(i,i);
            temp(i) = u(i) - (w/a)*(dot(A(i,1:i-1),temp(1:i-1)) + dot(A(i,i:n),u(i:n)) - f(i));
        end
        for i = fliplr(1:n)
            a = A(i,i);
            u(i) = temp(i) - (w/a)*(dot(A(i,i+1:n),u(i+1:n)) + dot(A(i,1:i),temp(1:i)) - f(i));
        end  
        err = max(abs(A*u-f));        
    end
    

end
