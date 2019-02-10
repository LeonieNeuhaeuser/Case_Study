
function [R P A_bar] =multigrid_operator(A)

% Reference: Prolongation &Restriction from NLA HW Q2
%                  Matrix coarse operator from NLA Lecture note

% n is the dimension of the  input vector which should have 
% The dimension of A mus be odd 

n = size(A,2);





e = ones(n,1);
P = spdiags([e/2 e e/2],-1:1,n,n);
P = P(:,2:2:end);

R = 1/2* P';

A_bar =  R*A*P;

end

