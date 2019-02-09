function [uk,errvect,iterations]=Conjugate_gradient(A,RHS,u0,tol)
%pk is the search direction
%rk and rk_plus_1 are the residuals
%Apk is A*pk

%initialize
rk=RHS-A*u0;
pk=rk;
Apk=A*pk;
norm_rk=norm(rk,2);
errvect(1)=norm_rk;
%loop
k=1;
uk=u0;
while norm_rk>tol
    alfa_k=norm_rk^2/((pk.')*Apk);
    uk=uk+alfa_k*pk;
    rk_plus_1=rk-alfa_k*Apk;
    
    norm_rk_plus_1=norm(rk_plus_1,2);
    
    beta_k_plus_1=norm_rk_plus_1^2/norm_rk^2;
    pk=rk_plus_1+beta_k_plus_1*pk; %in essence this is pk+1 
    Apk=A*pk;
    %but we can overwrite as we don;t need 2 lvls of pk at the same time
    k=k+1;
    norm_rk=norm_rk_plus_1;
    rk=rk_plus_1;
    
    errvect(k)=norm_rk;
end
iterations=k-1;
end