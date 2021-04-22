% % m, n, k
% % Manifold of m-by-n real matrices of fixed rank k 
% %       R is (m - r)x(m)
% % Where m = (ell+1)*(p_out + m_in)
% % Manifold: m = (ell+1)*q - r
% %           n = (ell+1)*q 
% %           k = r = (ell + 1)* m_in + ell * p_out
% k = r;
k_dim = (ell + 1)* m_in + ell * p_out;
m_dim = (ell+1)*q-k_dim;
n_dim = (ell+1)*q;
M_dim = (m_dim+n_dim-k)*k;
Delta0 = sqrt(M_dim) / 8 ;
Delta = Delta0;
% % Random vector in T_x M (this has to be very small)
% %%
% Rrand = rand(2,2) * diag(sort(rand(k, 1), 1, 'descend'))  * rand(2,12);
% Rrand2 = rand1() * diag(sort(rand(10, 1), 1, 'descend')) * rand2()';
% 
% Rrand * Rrand' 
% Rrand2 * Rrand2' 
% 
% %%
% rand1 = @() qr_unique(randn(2, 10));
% rand2 = @() qr_unique(randn(12, 10));

%%
eta = lincomb(RmyOpt, 1e-6, randomvec(RmyOpt));
%%
% Must be inside trust-region
%
while norm(eta(:)) > Delta
    eta = lincomb(x, sqrt(sqrt(eps)), eta);
end

%% STIEFEL
% R^nxp
% function M = stiefelfactory(n, p)
% X'*X = eye(p)
n_stiefel = size(Rini,1);
p_stiefel = size(Rini,2);

x = R'
R = x' 
%%
egrad = slra_mex_obj('grad', obj, RmyOpt)'  
%%
ehess = slra_mex_obj('mhess', obj, RmyOpt, eta);
rhess = ehess2rhess(RmyOpt', egrad, slra_mex_obj('mhess', obj, RmyOpt, eta)', eta')'

rhess2 = hessM(RmyOpt, eta, obj)*2

diff_rhess = rhess - rhess2;
norm(rhess(:))
norm(rhess2(:))
norm(diff_rhess(:))
(norm(diff_rhess(:))) / norm(rhess(:)) * 100

%%

(rhess + eye(size(rhess,2)))*egrad'
%%
function rhess = ehess2rhess(X, egrad, ehess, H)
    % Euclidean part
    XtG = multiprod(multitransp(X), egrad);
    symXtG = multisym(XtG);
    HsymXtG = multiprod(H, symXtG);
    rhess = projection(X, ehess - HsymXtG);
end

function Hh = hessM(U, H, obj)
%     fdeps = 1e-12;
    fdeps = 1;
    G0 = slra_mex_obj('grad', obj, U);
    G1 = slra_mex_obj('grad', obj, (U+fdeps*H));
    Hh = projection(U, (G1 - G0)/fdeps);
  end
% % Linear combination of tangent vectors
function v = lincomb(x, a1, d1, a2, d2)

    if nargin == 3
        v = a1*d1;
    elseif nargin == 5
        v = a1*d1 + a2*d2;
    else
        error('matrixlincomb takes either 3 or 5 inputs.');
    end

end

function U = randomvec(X)
%     U = projection(X, randn(n_stiefel, p_stiefel));
    U = projection(X, randn(2, 12));
    U = U / norm(U(:));
end
