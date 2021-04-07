image=zeros(1000,1000,3); %initialize
image(:,:,1)=1;   %Red (maximum value)
image(:,:,2)=rand(1000,1000);  %Green
image(:,:,3)=rand(1000,1000);  %Blue 
figure, imshow(image)

%%
clear
n = 100;
bnds = randn(n,2);
l = min( bnds, [], 2 );
u = max( bnds, [], 2 );


%% 
% Generate random problem data.
n = 1000;
A = randn(n);
A = .5*(A+A.');
 
% Create the problem structure.
manifold = spherefactory(n);
problem.M = manifold;
 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) -x'*(A*x);
problem.egrad = @(x) -2*A*x;      % notice the 'e' in 'egrad' for Euclidean
 
% Numerically check gradient consistency (optional).
checkgradient(problem);
 
% Solve.
[x, xcost, info, options] = trustregions(problem);
 
% Display some statistics.
figure;
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iteration number');
ylabel('Norm of the gradient of f');