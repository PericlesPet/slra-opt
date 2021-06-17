%% INITIATE SLRA
% Generally, series of operations:
% slra_sw_example -> opt_algos/varpro_evaluation_test -> opt_algos/augLag/alm/almTest -> opt_algos/myFmincon 

    % ADD PATH
addpath(genpath('..\..\..\Matlab\slra-slra-b1908bf'))
addpath(genpath('..\..\..\Matlab\slra-opt'))
clc
clear all, randn('seed', 0), rand('seed', 0), warning off

    % DESIGN PARAMETERS
m_in = 2;            % Inputs
p_out = 2;           % Outputs
ell = 2;             % l time-horizon / dynamics                                                                                                  
s_noise = 0.10;      % Noise Variation

    % GENERATE RANDOM SYSTEM
[n, q, T, sys0, u0, y0, w0, u, y, w] = generate_random_sys(m_in, p_out, ell, s_noise);

    % SLRA Solver Options
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
% profile on
% profiler_data = profile('info');
%   GET IDENT SOLUTION
% tic, sysh_ident = ident(w, m, ell, opt_oe); t_ident = toc;
tic, sysh_ident = ident_custom(w, m_in, ell, opt_oe); t_ident = toc;
%   GET KUNG REALIZATION SOLUTION
tic, sysh_kung = w2h2ss(w, m_in, n); t_kung = toc;


    % GET ACCURACIES
% Confirm sys0 is the ground truth
sys_comparison(u0, y0, r2ss(ss2r(sys0), m_in, ell), 0);
% How Close is initial system to noisy [u,y] ??
M_noise = sys_comparison(u, y, sys0, 0);
% How close is SLRA IDENTIFIED system to initial system
M_ident = sys_comparison(u0, y0, sysh_ident, 0, t_ident);
% How close is KUNG REALIZATION system to initial system
M_kung = sys_comparison(u0, y0, sysh_kung, 0, t_kung);




%% Perform up to "Data Generation" of slra_sw_example
tic
[w, s, r, opt, q, N, T] = sys2slra(w, m_in, ell, opt_mo);
toc

p = w2p(w);
Rh = ss2r(sys0); 
np = length(p);     

tic
[tts, p_new, p, r, s, w_new, Rini, phi ,psi, opt, ...
    th2R, C, s0, prob, pext, bfs, wtfdata] = ... 
    slra2slra_ext(p, s, r, opt);
toc

    % SLRA_MEX_OBJ
obj = slra_mex_obj('new', p, s, r);
if ~exist('info') || ~exist('info_ext')
    tic, [ph, info] = slra_mex_obj('optimize', obj, opt); t_slra = toc;
    tic, [ph_ext, info_ext, slraProbeData] = slra_ext(s2s(s, np), p, r, s.w, opt.Rini, s.phi, opt.psi, opt); t_slraext = toc;
    % probeData is just to check specific things
    toc
end

    % Set various R's
if exist('sysh_ident', 'var'), R_ident = ss2r(sysh_ident); end
if exist('sysh_kung', 'var'), R_kung = ss2r(sysh_kung); end
if exist('sys0', 'var'), R_true = ss2r(sys0); end
R_slramex = info.Rh;
R_slraext = info_ext.Rh;
Rin = Rini;
    % Set Helper Functions
sysAccuracy = @(R)compare(iddata(slradata.y0,slradata.u0),idss(r2ss(R,slradata.m_in,slradata.ell)));
p_R2X = @(p, R) [p(:); R(:)];
p2pext      = @(p) [0 ; p];
pext2hankel = @(pext) pext(wtfdata.tts+1);
x2hankel    = @(X) (wtfdata.s0 + pext2hankel(p2pext(X(1:np))));
x2R         = @(X) (reshape(X(np+1:end), size(Rini)));

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

