%% INITIATE SLRA
% Generally, series of operations:
% slraInit -> opt_algos/gd/gdMain
%          -> opt_algos/augLag/alm/almTest 
%          -> opt_algos/fminconTests 
%          -> opt_algos/panoc
%          -> opt_algos/fminlfbgs

    % ADD PATH
addpath(genpath('..\..\..\Matlab\slra-slra-b1908bf'))
addpath(genpath('..\..\..\Matlab\slra-opt'))
clc
clear all, randn('seed', 0), rand('seed', 0), warning off

%%
clc

    % DESIGN PARAMETERS
m_in = 2;            % Inputs
p_out = 2;           % Outputs
ell = 2;             % l time-horizon / dynamics                                                                                                  
s_noise = 0.1;      % Noise Variation
multiplier = 5;      % Multiplies base amount of samples to generate

M_ident = 0;
while M_ident <= 75
        % GENERATE RANDOM SYSTEM
    [n, q, T, sys0, u0, y0, w0, u, y, w] = generate_random_sys(m_in, p_out, ell, s_noise, multiplier);
        % SLRA Solver Options
    complexity = T*m_in*p_out;
    fprintf('Running for Complexity: %d\n\n', T*m_in*p_out);

    opt_oe.exct = 1:m_in;      % fixed inputs = output error identification
    use_c = 1;
    if use_c
        opt_oe.solver = 'c';
    else
        opt_oe.solver = 'm';
        opt_oe.method = 'reg';
    end
    % opt_MyOptimization
    opt_mo.exct = 1:m_in;      % fixed inputs = output error identification
    opt_mo.solver = 'm';
    opt_mo.method = 'reg';

        % IDENT
        % GET IDENT SOLUTION
    % tic, sysh_ident = ident(w, m, ell, opt_oe); t_ident = toc;
    tic, [sysh_ident, info_ident, wh_ident, xini_ident] = ident_custom(w, m_in, ell, opt_oe); t_ident = toc;

    % GET M IDENT
    getMident = 0;
    if getMident
        tic, [syshM_ident, infoM_ident, whM_ident, xiniM_ident] = ident_custom(w, m_in, ell, opt_mo); t_Mident = toc;
    end
    %   GET KUNG REALIZATION SOLUTION
    tic, sysh_kung = w2h2ss(w, m_in, n); t_kung = toc;


        % GET ACCURACIES
    % Confirm sys0 is the ground truth
    fprintf('Ground Truth (Check): \n');
    sys_comparison(u0, y0, r2ss(ss2r(sys0), m_in, ell), 0);
    % How Close is initial system to noisy [u,y] ??
    fprintf('Ground Truth + Noise: \n');
    M_noise = sys_comparison(u, y, sys0, 0);
    % How close is SLRA IDENTIFIED system to initial system
    fprintf('Ident system: \n');
    M_ident = sys_comparison(u0, y0, sysh_ident, 0, t_ident);
    % How close is KUNG REALIZATION system to initial system
    fprintf('Kung Realization: \n');
    M_kung = sys_comparison(u0, y0, sysh_kung, 0, t_kung);

    if M_ident <= 75
        clc
        fprintf('REPEAT\n');
    end
end

%
fprintf('First Part of SLRA data generation COMPLETE\n')

    % MIDDLE
    % Perform up to "Data Generation" of slra_sw_example
tic
[w, s, r, opt, q, N, T] = sys2slra(w, m_in, ell, opt_mo);
fprintf('sys2slra complete, t = %f\n', toc);

p = w2p(w);
Rh = ss2r(sys0); 
np = length(p);     

tic
[tts, p_new, p, r, s, w_new, Rini, phi ,psi, opt, ...
    th2R, C, s0, prob, pext, bfs, wtfdata] = ... 
    slra2slra_ext(p, s, r, opt);
fprintf('slra2slra_ext complete, t = %f\n', toc);

    % SLRA_MEX_OBJ
obj = slra_mex_obj('new', p, s, r);

if ~exist('info_mex') || ~exist('info_ext')
    tic, [ph, info_mex] = slra_mex_obj('optimize', obj, opt); t_slra = toc;
   
    %     tic, [ph_ext, info_ext, slraProbeData] = ...
    %     slra_ext(s2s(s, np), p, r, s.w, opt.Rini, s.phi, opt.psi, opt); t_slraext = toc;
    % probeData is just to check specific things
end

    % Set various R's
if exist('sysh_ident', 'var'), R_ident = ss2r(sysh_ident); end
if exist('sysh_kung', 'var'), R_kung = ss2r(sysh_kung); end
if exist('sys0', 'var'), R_true = ss2r(sys0); end
R_slramex = info_mex.Rh;    % R_slraext = info_ext.Rh;
Rin = Rini;

    % Set Main Structs
slradata.obj     = obj;
slradata.np      = length(p);
slradata.npExt   = length(wtfdata.p);
slradata.p       = p;
% slradata.ce     = problem.ce;
slradata.Rini    = Rini;
slradata.m_in    = m_in;
slradata.p_out   = p_out;
slradata.ell     = ell;
slradata.y0      = y0;
slradata.u0      = u0;
slradata.s       = s;
slradata.wtfdata = wtfdata;
slradata.f_opt   = slra_mex_obj('func', obj, R_ident);

    % Set Helper Functions
sysAccuracy = @(R)compare(iddata(slradata.y0,slradata.u0),idss(r2ss(R,slradata.m_in,slradata.ell)));
p_R2X = @(p, R) [p(:); R(:)];
p2pext      = @(p) [0 ; p];
pext2hankel = @(pext) pext(wtfdata.tts+1);
x2hankel    = @(X) (wtfdata.s0 + pext2hankel(p2pext(X(1:np))));
x2R         = @(X) (reshape(X(np+1:end), size(Rini)));




    % Set Various Other

for almInit = 1
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
        'niter', 10 , ... 
        'miter', 100 ...
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
            'rho', checkdata.rho(end), ...        
            'niter', 20 , ... 
            'miter', 100 ...
        );
        x0 = checkdata.x(:,end);
        lambda0 = checkdata.lambda(:,end);
        checkdata_temp = checkdata;
    end
     

    
end

f_slra = slradata.f_opt;
f_opt = slra_mex_obj('func', obj, R_true);

constraints.ce          = ce;
constraints.dce         = dce;
constraints.ce_rri2     = @(x) stiefConstraint(reshape(x,size(Rini)), 'dist');
constraints.dce_rri2    = @(x)FDGradient(constraints.ce_rri2,x,gstep);

[~, M_slra] = sysAccuracy(R_ident);
slradata.M0 = M_slra;
dimensions  = size(Rini);
dimension   = dimensions(1) * dimensions(2);
LM_obj_reg  = @(R) SLRA_FuncAndGrad_vals(obj, opt.g, dimensions, R, 'reg');
% get starting value gamma
beta_safety_value=0.05;
df  = @(x) reshape(slra_mex_obj('grad', obj, reshape(x, dimensions)), dimension , 1);
lipschitz_constant = estimate_lipschitz( df,Rini(:));
gamma0 = (1-beta_safety_value)/lipschitz_constant;
gamma  = gamma0;

    % Get ident Accuracies
tic
if ~exist('Ms_ident')
    for i = 1:length(info_ident.iterinfo)
       R_curr = info_ident.RhK(:,:,i);
       [~, M_curr] = sysAccuracy(R_curr);
       Ms_ident(:, i) = M_curr;
    end
end
t_mIdent = toc;

plotM_ident = 0;
if plotM_ident
    figure
    plot( info_ident.iterinfo(1,:), mean(Ms_ident))
    
    figure
    for i = 1:size(Ms_ident,1)
        plot( info_ident.iterinfo(1,:), Ms_ident(i,:))
        hold on
    end
end
fprintf('\nIdent Accuracies Obtained, t = %f\n', t_mIdent)
fprintf('slra has been initiated, proceed to optimization algos\n')

isCloseAll   = 1;

%% STEP RESPONSES
% Gather All Systems Into a Struct
clear systems;
systems.sys0        = sys0;
systems.sysh_ident  = sysh_ident;
systems.sysh_kung   = sysh_kung;

% Plot Step Responses
systems_step_responses(systems)



%% DELETE OBJ
slra_mex_obj('delete', obj);

