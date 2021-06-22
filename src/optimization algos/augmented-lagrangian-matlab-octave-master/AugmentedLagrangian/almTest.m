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

    % Get weights and indices
    alm_weights         = s.w;
    non_inf_indices     = ~isinf(alm_weights);
    inf_indices         = isinf(alm_weights);

    % Get x0, Transform to reduced_x0
    x   = [ph_ini(:) ; Rini(:)];
    x0  = x;
	
    np_w = length(wtfdata.p);
    % Helper lambda functions
    p2pext      = @(X) [0 ; X];
    pext2hankel = @(X) X(wtfdata.tts+1);
    x2hankel    = @(X) (wtfdata.s0 + pext2hankel(p2pext(X(1:np_w))));
    x2R         = @(X) (reshape(X(np_w+1:end), size(Rini)));
       
    % Function and Gradient 
    %     alm_weights = [1e6*ones(sum(isinf(s.w)), 1); ones(sum(~isinf(s.w)), 1)];
    f = @(X) 1/2*norm(alm_weights(~isinf(alm_weights)).* (X(1:np_w)-wtfdata.p(1:np_w)))^2 ; 
    df = @(X) [alm_weights(~isinf(alm_weights)).* (X(1:np_w)-wtfdata.p(1:np_w)) ; zeros(R_n * R_m, 1)] ;

    
    % Constraint equality with norm
    gstep = 1/1e9;
    % norm( R*H(X) ) = 0  Constraint
    ce_rh0 = @(X) norm((reshape (... 
            x2R(X) * x2hankel(X), ... 
            size(Rini, 1) * (T-ell), 1 ...
            )));
    % norm( R*R - I ) = 0 Constraint
    ce_rri = @(X) stiefConstraint(x2R(X), 'dist');

    % (Forward) Gradients of Constraints
    dce_rri = @(X) FDGradient(ce_rri, X, gstep);
    dce_rh0 = @(X) FDGradient(ce_rh0, X, gstep);
  
    % All Constraints in 1 
    ce = @(X)([ce_rh0(X), ce_rri(X)]);
    dce = @(X) ([dce_rh0(X) dce_rri(X)]);
    % 'f', @(X) 1/2*norm(X(1:np)-p)^2, ...
    % 'df', @(X) [(X(1:np)-p) ; zeros(R_n * R_m, 1)] , ...
    problem = struct( ...
        'f', f, ...
        'df', df, ...
        'ce', ce, ...
        'dce', dce ... 
        );
        % 'ci', @(X) (0), ...
        % 'dci', @(X) 0 ...

    lambda0 = [];

    options = struct( ...
        'gma', 5 ,...
        'niter', 1 , ... 
        'miter', 1 ...
    );

    %     slradata.obj    = obj;
    %     slradata.np     = np;
    %     slradata.p      = p;
    slradata.ce     = problem.ce;
    %     slradata.Rini   = Rini;
    %     slradata.m_in   = m_in;
    %     slradata.ell    = ell;
    %     slradata.y0     = y0;
    %     slradata.u0     = u0;
    %     slradata.s      = s;

    continue_from_earlier = 0;
    if continue_from_earlier == 1
        options = struct( ...
            'gma', 5 ,...
            'rho', almData.rho(end), ...        
            'niter', 20 , ... 
            'miter', 100 ...
        );
        x0 = almData.x(:,end);
        lambda0 = almData.lambda(:,end);
        checkdata_temp = almData;
    end
     

    
 end

tic
[x, fval, lambda, kkt, almData] = ...
    almSolve(problem, x0, lambda0, options, slradata);
toc

doVisual = 0;
if doVisual == 1
    almVisualize
end

% % % Get final x
% % % x_final(inf_indices,1)        = p(inf_indices);
% % % x_final(non_inf_indices,1)    = x(1:np);
% % % 
% % % R_alm_vec = x(np_w+1:end);
% % % R_alm = reshape(R_alm_vec, R_n, R_m);
% % % ph_alm = x(1:np_w);
% % % 
% % % f_temp = slra_mex_obj('func', obj, R_alm);
% % % f_temp = slra_mex_obj('func', obj, Rini);
% % % [~, M_alm] = accuracy_r(R_alm);
% % % [~, M_ini] = accuracy_r(Rini);
% % % [~, M_slraOpt] = accuracy_r(info1.Rh);
% % % mean([M_alm M_ini M_slraOpt]);
% % % 
