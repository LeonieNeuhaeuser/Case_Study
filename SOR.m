function [u, err, errvec] = SOR(w, A, f, u0, tol)

    n = length(u0);
    L = tril(A, -1);
    R = triu(A, 1);
    D = spdiags(diag(A), 0, n, n);
    errvec =[];
    u = u0;
    err = max(abs(A*u - f));
    
    while err > tol
        u = (D+w*L)\(((1-w)*D-w*R)*u + w*f);
        err = max(abs(A*u - f));
        %disp(err);
        errvec =[errvec err];
    end
end