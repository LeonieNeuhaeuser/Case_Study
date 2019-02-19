%% Error plots 
clear, close all

% Example problems 
f = @(x,y) 13*pi^2*sin(2*pi*x).*sin(3*pi*y);
g = @(x,y) -(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5; 
sol1 = @(x,y) sin(2*pi*x).*sin(3*pi*y);
sol2 = @(x,y) (x-1).^5.*x.^2.*y.*(y-1);  

% Parameters 
xmin = 0; xmax = 1; ymin = 0; ymax = 1;
tol = 1e-5;
presteps = 5;
poststeps = 5; 
levels = 5; 
nmin = 1000;

%but are preferred for SOR and SSOR 
N = 512; 
M = 512; 

% Initialise 
dx = 1/N; dy = 1/M;
x = xmin:dx:xmax;
y = ymin:dy:ymax;
[X,Y] = meshgrid(x,y); 
u1 = zeros(M+1,N+1); 
u2 = zeros(M+1,N+1);

% Descretising the system
A = create_2d_finite_diff_A(xmin,xmax,ymin,ymax,N,M);
F = makeF(X,Y,f,N,M);%generate the RHS
G = makeF(X,Y,g,N,M);

% Solving the system and reshaping
u0 = zeros(length(F),1);

%start tic 2
% Pre & Post smoothing step set to be 5 repetively
disp('Started multigrid')
tic
[y1_MG, errvect_1] = Multigrid_Vcycle(A,F,u0,tol,presteps,poststeps,nmin,1,levels);
toc
%endtic2
%starttic 3
tic
[y2_MG, errvect_2] = Multigrid_Vcycle(A,G,u0,tol,presteps,poststeps,nmin,1,levels);
toc
%endtic3

y_1 = reshape(y1_MG, M-1, N-1); u1(2:M,2:N) = y_1;
y_2 = reshape(y2_MG, M-1, N-1); u2(2:M,2:N) = y_2;

% Plotting the errors
figure
semilogy(errvect_1)
title('Test problem 1')
ylabel('error') 
xlabel('iterations') 

figure
semilogy(errvect_2)
title('Test problem 2')
ylabel('error') 
xlabel('iterations') 

%%
xx=linspace(0,1,N+1);
yy=linspace(0,1,M+1);
[X,Y]=meshgrid(xx,yy);
figure
surf(X,Y,u1)
title("Numerical Solution 1 - Multigrid");
xlabel("X axis")
ylabel("Y axis")
zlabel("Z axis")