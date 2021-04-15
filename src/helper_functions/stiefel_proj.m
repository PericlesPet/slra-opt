% STIEFEL PROJECTION
% According to: 
% Optimization Algorithms on Matrix Manifolds
% By Authors : P.-A. Absil , R. Mahony , R. Sepulchre
% Page ~81

% P_x(Z) = Z - X * sym(X'Z)
% Where sym(M) denotes the symmetric part of M,
%     sym(M) = 1/2 * (M + M')

% For Stiefel Projected Gradient nabla_f(x) 
% of a function f_hat(x) subject to x on a Stiefel Manifold,
% nabla_f(x) = nabla_f_hat(x) - x*symmetric(x' * nabla_f_hat(x))
% Therefore:
% nabla_f_hat(x) --> Z
% x              --> X
function R_stiefel = stiefel_proj(X, Z)

R_stiefel = Z - X*symmetric(X'*Z);

end