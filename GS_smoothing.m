function u = GS_smoothing(A,f, u0, steps, w)
    
% Using relaxed Gauss-Seidel smoothing
    n = length(u0);
    L = tril(A, -1);
    R = triu(A, 1);
    D = spdiags(diag(A), 0, n, n);
    u = u0;
    
    for i = 1:steps
        u = (D+w*L)\(((1-w)*D-w*R)*u + w*f);
    end
    
end 