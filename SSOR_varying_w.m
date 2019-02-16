
%% Error plots 
clear, close all
omegavector = linspace(0,2,11);
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
    [y_1_SSOR, err_1_SSOR, errvec_1_SSOR] = SSOR2(w,A,F,u0,tol);
    semilogy([1:numel(errvec_1_SSOR)],errvec_1_SSOR, 'DisplayName',sprintf('%i',w));
    legend('show')
    hold on
end
title('Convergence plots for different omega - Sym-SOR- Test 1')

figure
for i=1:numel(omegavector)
    w=omegavector(i);
    [y_2_SSOR,err_2_SSOR, errvec_2_SSOR] = SSOR2(w,A,G,u0,tol);
    semilogy([1:numel(errvec_2_SSOR)],errvec_2_SSOR, 'DisplayName',sprintf('%i',w));
    legend('show')
    hold on
end
title('Convergence plots for different omega - Sym-SOR- Test 2')

%% Convergence rates as function of w 

omegavector = linspace(0.01,1.99,200); 
k_1 = zeros(1,length(omegavector)); k_2 = k_1; 

% Solving the system and reshaping
u0 = zeros(length(F),1);
for i = 1:numel(omegavector)
    w = omegavector(i);
    [y_1_SSOR, err_1_SSOR, errvec_1_SSOR] = SSOR2(w,A,F,u0,tol);
    k_1(i) = length(errvec_1_SSOR); 
end
figure
semilogy(omegavector, k_1)
xlabel('$\omega$','Interpreter','latex')
ylabel('Iterations', 'Interpreter', 'latex')

for i = 1:numel(omegavector)
    w = omegavector(i);
    [y_2_SSOR,err_2_SSOR, errvec_2_SSOR] = SSOR2(w,A,G,u0,tol);
    k_2(i) = length(errvec_2_SSOR); 
end
figure
semilogy(omegavector, k_2)
xlabel('$\omega$','Interpreter','latex')
ylabel('Iterations', 'Interpreter', 'latex')
