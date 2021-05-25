% Simple problems from wikipedia and elsewhere for testing the Augmented Lagrangian optimization.
%
% Implementation tips:
% ** Nonlinear programming problems tend to converge more reliably when the initial point, x0, is
%   feasible (feasibility of the initial point, however, is not required by the algorithm).
% ** Another strategy that may be used to improve convergence on nonlinear problems it to initialize
%   the Lagrange multipliers using the Moore-Penrose pseudoinverse:
%       lambda0 = pinv([dce(x0), dci(x0)])*df(x0)
%   which is the default initialization when lambda0 is not provided by the user (but with slack
%   variables included).
% ** Nonlinear problems may also require tuning of the increment parameter, gma, to improve
%   convergence.

p_alm = 0;

gma = 10;

 for alm_problem_init = 1
    % minimize f(x), subject to ce(x) = 0, ci(x) >= 0 
    % Minimize   : f(X) = f(ph, R) = 1/2 * norm(ph - p)^2, 
    % Subject To : R * Hankel(ph) = 0  (... later ... , R' * R = I)

    % PARAMETERS
    m_t  = (ell+1) * (m_in + p_out);
    n_t  = (T-ell);
    np = length(p);
    phi  = eye(m_t);
    pext = [0; p];
    tts  = s2s(s,np);
    vec_tts = tts(:); NP = 1:np;
    bfs = vec_tts(:, ones(1, np)) == NP(ones(m_t * n_t, 1), :);
    s0_2   = zeros(m_t, n_t);

    R_n = size(Rini, 1);
    R_m = size(Rini, 2);

    mslra_handle = @(R) Mslra_ext(R, wtfdata.tts, wtfdata.p, [], ... 
        wtfdata.bfs, wtfdata.phi, wtfdata.s0);

    % Get initial P_hat, x0
    [~, ph_ini] = mslra_handle(Rini); 
    ph_ini = [p(isinf(s.w)); ph_ini];
    X = [ph_ini(:) ; Rini(:)];
    x0 = X;

    % Constraint equality with norm
    gstep = 1/1e9;
    ce_n = @(X) norm((reshape (... 
            reshape(X(np+1:end), size(Rini)) * X(tts), ... 
            size(Rini, 1) * (T-ell), 1 ...
            )));
    dce_n = @(X) FDGradient(ce_n, X, gstep);

%         'f', @(X) 1/2*norm(X(1:np)-p)^2, ...
%         'df', @(X) [(X(1:np)-p) ; zeros(R_n * R_m, 1)] , ...
    problem = struct( ...
        'f', @(X) 1/2*norm(X(np/2+1:np)-p(np/2+1:np))^2, ...
        'df', @(X) [zeros(np/2, 1) ; (X(np/2+1:np)-p(np/2+1:np)) ; zeros(R_n * R_m, 1)] , ...
        'ce', ce_n, ...
        'dce', dce_n ... 
        );
% 		'ci', @(X) (0), ...
% 		'dci', @(X) 0 ...

    lambda0 = [];


    slradata.obj    = obj;
    slradata.np     = np;
    slradata.p      = p;
    slradata.ce     = problem.ce;
    slradata.Rini   = Rini;
    slradata.m_in   = m_in;
    slradata.ell    = ell;
    slradata.y0     = y0;
    slradata.u0     = u0;
    slradata.s      = s;
 end

%  'gma', gma ,...
options = struct( ...
    'gma', 15 ,...
    'niter', 8 , ... 
    'miter', 30 ...
);


tic
[x, fval, lambda, kkt, checkdata] = ...
    almSolve(problem, x0, lambda0, options, slradata);
toc
%%
if ~exist('info1')
    tic
    [ph1, info1] = slra(p, s, r, opt)
    toc
    accuracy_r = @(R)compare(iddata(slradata.y0, slradata.u0), idss(r2ss(R, slradata.m_in, slradata.ell))); 
end

almVisualize


%%
R_alm_vec = x(np+1:end);
R_alm = reshape(R_alm_vec, R_n, R_m);
ph_alm = x(1:np);

R_alm = reshape(x(np+1:end), 2, 12);
f_temp = slra_mex_obj('func', obj, R_alm)
f_temp = slra_mex_obj('func', obj, Rini)
[~, M_alm] = accuracy_r(R_alm);
[~, M_ini] = accuracy_r(Rini);
[~, M_slraOpt] = accuracy_r(info1.Rh);
mean([M_alm M_ini M_slraOpt])
% sys_comparison(u0, y0, r2ss(Rini, m_in, ell))

