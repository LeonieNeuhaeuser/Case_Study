%conjugate gradients error plots for varying mesh size
%% mesh evaluation of Conjugate_gradient
clear, close all

% Example problems 
f = @(x,y) 13*pi^2*sin(2*pi*x).*sin(3*pi*y);
g = @(x,y) -(x-1).^3.*(42*x.^2-24*x+2).*y.*(y-1)-2*x.^2.*(x-1).^5; 
sol1 = @(x,y) sin(2*pi*x).*sin(3*pi*y);
sol2 = @(x,y) (x-1).^5.*x.^2.*y.*(y-1);  

% Parameters 
xmin = 0; xmax = 1; ymin = 0; ymax = 1;
tol = 1e-3;
%w=1
iter_vec1=[];
iter_vec2=[];

for N= [8,16,32,64,128,256,512]

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

        [y_1_CG, errvec_1_CG,iter1] = Conjugate_gradient(A,F,u0,tol);
        y_1_CG = reshape(y_1_CG, N-1, N-1); u1(2:N,2:N) = y_1_CG;
        iter_vec1 = [iter_vec1 iter1];

    % Example Problem 2

        [y_2_CG,errvec_2_CG,iter2] = Conjugate_gradient(A,G,u0,tol); 
        y_2_CG = reshape(y_2_CG, N-1, N-1); u2(2:N,2:N) = y_2_CG;
        iter_vec2 = [iter_vec2 iter2];

end

h=figure
    semilogy([8,16,32,64,128,256,512],iter_vec1)
    title('Test problem 1')
    xlabel('Mesh Size')
    ylabel('Number of Iterations')

g=figure
    semilogy([8,16,32,64,128,256,512],iter_vec2)
    title('Test problem 2') 
    xlabel('Mesh Size')
    ylabel('Number of Iterations')
