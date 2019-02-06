function A=create_2d_finite_diff_A(a,b,c,d,N,M)
%N-number of pts in x dir
%M number of pts in Y dir
%ab is the x interval
%cd is the y interval
dx=abs(b-a)/N;
dy=abs(c-d)/M;
e1=  repmat(-1/dy^2,(M-1)*(N-1),1);
e2=repmat(2*(dx^(-2)+dy^(-2)),(M-1)*(N-1),1);
e3=repmat(-1/dx^2,(M-1)*(N-1),1);
A = spdiags([e3 e1 e2 e1 e3 ], [-M+1,-1,0,1,M-1],(M-1)*(N-1), (M-1)*(N-1));
for k=1:N-2
    A(k*(M-1),k*(M-1)+1)=0;
    A(k*(M-1)+1,k*(M-1))=0;
end
end