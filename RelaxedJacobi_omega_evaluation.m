%% Omega evaluation of Relaxed Jacobi 
clear, close all

% Example problems 
f = @(x,y) 13*pi^2*sin(2*pi*x).*sin(3*pi*y);
g = @(x,y) -(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5; 
sol1 = @(x,y) sin(2*pi*x).*sin(3*pi*y);
sol2 = @(x,y) (x-1).^5.*x.^2.*y.*(y-1);  

% Parameters 
xmin = 0; xmax = 1; ymin = 0; ymax = 1;
tol = 1e-3;
N = 32; 
M = 32; 

% Initialise 
dx = 1/N; dy = 1/M;
x = xmin:dx:xmax;
y = ymin:dy:ymax;
[X,Y] = meshgrid(x,y); 
u1 = zeros(M+1,N+1); 
u2 = zeros(M+1,N+1);

% Descretising the system
A = create_2d_finite_diff_A(xmin,xmax,ymin,ymax,N,M);
F = makeF(X,Y,f,N,M);
G = makeF(X,Y,g,N,M);

% Solving the system and reshaping
u0 = zeros(length(F),1);
<<<<<<< HEAD
parameter = [0.1,0.5,1];
lambda = cos(pi*dx);
h = figure
=======
h = figure;
>>>>>>> bf97ab8f30666d75905dae39610c42618948c39a
% Example Problem 1
for w=parameter
    [y_1_RJ, err_1_RJ, errvec_1_RJ] = RelaxedJacobi(w,A,F,u0,tol); y_1_RJ = reshape(y_1_RJ, M-1, N-1); u1(2:M,2:N) = y_1_RJ;
    semilogy(errvec_1_RJ, 'Displayname', sprintf('%i',w));
    hold on
    len = length(errvec_1_RJ);
    curr_lam = (1-w)+w*lambda;
    x = linspace(0,len,len)
    semilogy(x,norm(u1(1:M+1,1:N+1))*curr_lam.^x)
    title('Test problem 1')
    legend('show');
end
hold off

g=figure;
% Example Problem 2
for w=parameter
    [y_2_RJ,err_2_RJ, errvec_2_RJ] = RelaxedJacobi(w,A,G,u0,tol); y_2_RJ = reshape(y_2_RJ, M-1, N-1); u2(2:M,2:N) = y_2_RJ;
    semilogy(errvec_2_RJ, 'Displayname', sprintf('%i',w))
    title('Test problem 2')
    hold on
    legend('show');
end
