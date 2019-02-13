%% mesh evaluation of Relaxed Jacobi 
clear, close all

% Example problems 
f = @(x,y) 13*pi^2*sin(2*pi*x).*sin(3*pi*y);
g = @(x,y) -(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5; 
sol1 = @(x,y) sin(2*pi*x).*sin(3*pi*y);
sol2 = @(x,y) (x-1).^5.*x.^2.*y.*(y-1);  

% Parameters 
xmin = 0; xmax = 1; ymin = 0; ymax = 1;
tol = 1e-3;
w = 1;
err_vec1=[];
err_vec2=[];
bound = [];


iterations = [8,16,32,64];

for N = iterations

    % Initialise 
    dx = 1/N; dy = dx;
    x = xmin:dx:xmax;
    y = ymin:dy:ymax;
    [X,Y] = meshgrid(x,y); 
    u1 = zeros(N+1,N+1); 
    u2 = zeros(N+1,N+1);

    % Descretising the system
    A = create_2d_finite_diff_A(xmin,xmax,ymin,ymax,N,N);
    F = makeF(X,Y,f,N,N);
    G = makeF(X,Y,g,N,N);

    % Solving the system and reshaping
    u0 = zeros(length(F),1);
    
    % Example Problem 1
    [y_1_RJ, err_1_RJ, errvec_1_RJ] = RelaxedJacobi(w,A,F,u0,tol); y_1_RJ = reshape(y_1_RJ, N-1, N-1); u1(2:N,2:N) = y_1_RJ;
    err_vec1 = [err_vec1 err_1_RJ];

    % Example Problem 2
    [y_2_RJ,err_2_RJ, errvec_2_RJ] = RelaxedJacobi(w,A,G,u0,tol); y_2_RJ = reshape(y_2_RJ, N-1, N-1); u2(2:N,2:N) = y_2_RJ;
    err_vec2 = [err_vec2 err_2_RJ];
    
    curr_bound = log(tol)/log(cos(pi*dx));
    bound = [bound curr_bound];
end

h=figure;
    plot(iterations,err_vec1)
    hold on
    plot(iterations,bound)
    legend("actual performance","upper bound", 'Location','North')
    title('Test problem 1')
    xlabel("number of gridpoints (N=M)")
    ylabel("iterations")

g=figure;
    plot(iterations,err_vec2)
    hold on
    plot(iterations,bound)
    legend("actual performance","upper bound", 'Location', 'North')
    title('Test problem 2')
    xlabel("number of gridpoints (N=M)")
    ylabel("iterations")
