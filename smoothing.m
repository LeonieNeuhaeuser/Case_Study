function u = smoothing(A,f,  u_i, steps,w)

% Using relaxed Jacobi for smoothing 


    n = length(u_i);
    I = speye(n);
    L = tril(A,-1); 
    R = triu(A,1); 
    D_inv = spdiags(1./diag(A),0,n,n);
    iterations = 0;
    
    u = u_i; 

    
    for i = 1:steps
        u = ((1-w).* I - w .* D_inv *(L+R))*u + w.*D_inv*f;
    end


end 