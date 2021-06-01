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
    [~, ph_ini, mslracheck] = mslra_handle(Rini); 

    % Get weights and indices
    alm_weights         = s.w;
    non_inf_indices     = ~isinf(alm_weights);
    inf_indices         = isinf(alm_weights);

    % Get x0, Transform to reduced_x0
    x   = [ph_ini(:) ; Rini(:)];
    x0  = x;
	
    np = length(wtfdata.p);
    % Helper lambda functions
    p2pext      = @(X) [0 ; X];
    pext2hankel = @(X) X(wtfdata.tts+1);
    x2hankel    = @(X) (wtfdata.s0 + pext2hankel(p2pext(X(1:np))));
    x2R         = @(X) (reshape(X(np+1:end), size(Rini)));
    
    % Function and Gradient 
%     alm_weights = [1e6*ones(sum(isinf(s.w)), 1); ones(sum(~isinf(s.w)), 1)];
    f = @(X) 1/2*norm(alm_weights(~isinf(alm_weights)).* (X(1:np)-wtfdata.p(1:np)))^2 ; 
    df = @(X) [alm_weights(~isinf(alm_weights)).* (X(1:np)-wtfdata.p(1:np)) ; zeros(R_n * R_m, 1)];

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
        'niter', 10 , ... 
        'miter', 100 ...
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
    
    continue_from_earlier = 0;
    if continue_from_earlier == 1
        options = struct( ...
            'gma', 5 ,...
            'rho', checkdata.rho(end), ...        
            'niter', 20 , ... 
            'miter', 100 ...
        );
        x0 = checkdata.x(:,end);
        lambda0 = checkdata.lambda(:,end);
        checkdata_temp = checkdata;
    end
     
    
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

% Get final x
x_final(inf_indices,1)        = p(inf_indices);
x_final(non_inf_indices,1)    = x(1:np);

if ~exist('info1')
    tic
    [ph1, info1] = slra(p, s, r, opt);
    toc
    accuracy_r = @(R)compare(iddata(slradata.y0, slradata.u0), idss(r2ss(R, slradata.m_in, slradata.ell))); 
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
end

test_alm_gd_visualize
