function [u,err, err_vec] = RelaxedJacobi(w,A,f,u0,tol)

    n = length(u0)
    I = eye(n);
    L = tril(A,-1); 
    R = triu(A,1); 
    D_inv = spdiags(1./diag(A),0,n,n); 
    size(A)
    
    
    u = u0; 
    err = max(abs(A*u-f))
    err_vec = []
    
    while err > tol 
        u = ((1-w).* I - w .* D_inv *(L+R))*u + w.*D_inv*f;
        err = max(abs(A*u-f));
        disp(err)
        err_vec = [err_vec,err]
    end

end

