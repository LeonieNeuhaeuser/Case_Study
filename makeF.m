function F = makeF(X,Y,f,N,M)

    F = f(X(2:N,2:M),Y(2:N,2:M)); 
    F = F(:); 

end