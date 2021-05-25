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

%         'gma', gma ,...
     options = struct( ...
        'gma', 15 ,...
        'niter', 8 , ... 
        'miter', 30 ...
    );

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
 
 
tic
 [x, fval, lambda, kkt, checkdata] = ...
     almSolve(problem, x0, lambda0, options, slradata);
toc

    % 
    % INNER ITERATION = 30, norm(x-x0) = 1.833717
    %    f = 68.7502,        f0 = 5.0072
    %    L = 955710.7884,    L0 = 0.2759
    %    CE = 4.3720,        CE0 = 0.0000
    %    DP = 1.2228,        DP0 = 0.7427
    % Elapsed time is 101.488286 seconds.

    
%%
tic
[ph1, info1] = slra(p, s, r, opt)
toc
    % STATS: 
    % iterations: 91
    %         funcCount: 2400
    %          stepsize: 0.0012
    %      lssteplength: 1
    %     firstorderopt: 0.3665
    %         algorithm: 'quasi-newton'
    %           message: 'Solver stopped prematurely.??fminunc stopped because it exceeded the function evaluation limit,?options.MaxFunctionEvaluations = 2400 (the default value).'
    %              fmin: 4.5845
    %                Rh: [2×12 double]
    %                Th: [24×94 double]
    %                 F: [1×94 double]
    % 
    % Elapsed time is 37.140498 seconds.

    
accuracy_r = @(R)compare(iddata(slradata.y0, slradata.u0), idss(r2ss(R, slradata.m_in, slradata.ell))); 
%%
test_alm_gd_visualize

 %%
 [logdata, data_opt, minf_log, x] = ...
    myALM_GDesc(f, df, x0, options)
