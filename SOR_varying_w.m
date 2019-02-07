
%% Error plots 
clear, close all
omegavector= linspace(0,2,11)
omegavector(1)=[];
omegavector(10)=[];
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
for i=1:numel(omegavector)
    w=omegavector(i);
    [y_1_SOR, err_1_SOR, errvec_1_SOR] = SOR(w,A,F,u0,tol);
    semilogy([1:numel(errvec_1_SOR)],errvec_1_SOR, 'DisplayName',sprintf('%i',w));
    legend('show')
    hold on
end
title('Convergence plots for different omega - SOR- Test 1')

figure
for i=1:numel(omegavector)
    w=omegavector(i);
    [y_2_SOR,err_2_SOR, errvec_2_SOR] = SOR(w,A,G,u0,tol);
    semilogy([1:numel(errvec_2_SOR)],errvec_2_SOR, 'DisplayName',sprintf('%i',w));
    legend('show')
    hold on
end
title('Convergence plots for different omega - SOR- Test 2')
