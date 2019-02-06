function [u,err] = SSOR2(w,A,f,u0,tol)

    n = length(u0);
    L = tril(A,-1); 
    R = triu(A,1); 
    D = spdiags(diag(A),0,n,n); 
    
    u = u0; 
    err = max(abs(A*u-f));
    
    while err > tol 
        temp = (D+w*L)\(((1-w)*D-w*R)*u + w*f);
        u = (D+w*R)\(((1-w)*D-w*L)*temp + w*f);
        err = max(abs(A*u-f));
    end

end
