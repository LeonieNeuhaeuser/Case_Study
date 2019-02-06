%% Error plots 
clear, close all

% Example problems 
f = @(x,y) 13*pi^2*sin(2*pi*x).*sin(3*pi*y);
g = @(x,y) -(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5; 
sol1 = @(x,y) sin(2*pi*x).*sin(3*pi*y);
sol2 = @(x,y) (x-1).^5.*x.^2.*y.*(y-1);  

% Parameters 
xmin = 0; xmax = 1; ymin = 0; ymax = 1;
tol = 1e-3;
w = 0.5; %values larger than 1 won't work for relaxed Jacobi, but are preferred for SOR and SSOR 
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
[y_1_RJ, err_1_RJ, errvec_1_RJ] = RelaxedJacobi(w,A,F,u0,tol); y_1_RJ = reshape(y_1_RJ, M-1, N-1); u1(2:M,2:N) = y_1_RJ;
[y_2_RJ,err_2_RJ, errvec_2_RJ] = RelaxedJacobi(w,A,G,u0,tol); y_2_RJ = reshape(y_2_RJ, M-1, N-1); u2(2:M,2:N) = y_2_RJ;
[y_1_GS, err_1_GS, errvec_1_GS] = GS2(A,F,u0,tol); y_1_GS = reshape(y_1_GS, M-1, N-1); u1(2:M,2:N) = y_1_GS;
[y_2_GS, err_2_GS, errvec_2_GS] = GS2(A,G,u0,tol); y_2_GS = reshape(y_2_GS, M-1, N-1); u2(2:M,2:N) = y_2_GS;
[y_1_SOR, err_1_SOR, errvec_1_SOR] = SOR(w,A,F,u0,tol); y_1_SOR = reshape(y_1_SOR, M-1, N-1); u1(2:M,2:N) = y_1_SOR;
[y_2_SOR,err_2_SOR, errvec_2_SOR] = SOR(w,A,G,u0,tol); y_2_SOR = reshape(y_2_SOR, M-1, N-1); u2(2:M,2:N) = y_2_SOR;
[y_1_SSOR, err_1_SSOR, errvec_1_SSOR] = SSOR2(w,A,F,u0,tol); y_1_SSOR = reshape(y_1_SSOR, M-1, N-1); u1(2:M,2:N) = y_1_SSOR;
[y_2_SSOR,err_2_SSOR, errvec_2_SSOR] = SSOR2(w,A,G,u0,tol); y_2_SSOR = reshape(y_2_SSOR, M-1, N-1); u2(2:M,2:N) = y_2_SSOR;

% Plotting the errors
figure
semilogy(errvec_1_RJ)
title('Test problem 1')
hold on
semilogy(errvec_1_GS)
semilogy(errvec_1_SOR)
semilogy(errvec_1_SSOR)
hold off
legend('RJ','GS','SOR','SSOR')
ylabel('error') 
xlabel('iterations') 

figure
semilogy(errvec_2_RJ)
title('Test problem 2')
hold on
semilogy(errvec_2_GS)
semilogy(errvec_2_SOR)
semilogy(errvec_2_SSOR)
hold off
legend('RJ','GS','SOR','SSOR')
ylabel('error') 
xlabel('iterations') 




