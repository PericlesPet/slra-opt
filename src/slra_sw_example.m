%% simulation parameters
% addpath(genpath('helper_functions'))
addpath(genpath('..\..\..\Matlab'))
clc
clear all, randn('seed', 0), rand('seed', 0), warning off
%%
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

%% compare ident
% profile on
% profiler_data = profile('info');
%   GET IDENT SOLUTION
% tic, sysh_ident = ident(w, m, ell, opt_oe); t_ident = toc;
tic, sysh_ident = ident_custom(w, m_in, ell, opt_oe); t_ident = toc;
%   GET KUNG REALIZATION SOLUTION
tic, sysh_kung = w2h2ss(w, m_in, n); t_kung = toc;



%%
% How Close is initial system to noisy [u,y] ??
sys_comparison(u, y, sys0)                  
% How close is SLRA IDENTIFIED system to initial system
sys_comparison(u0, y0, sysh_ident, t_ident) 
% How close is KUNG REALIZATION system to initial system
sys_comparison(u0, y0, sysh_kung, t_kung)   

%% STEP RESPONSES
% Gather All Systems Into a Struct
systems.sys0 = sys0;
systems.sys1 = sysh_ident;
systems.sys2 = sysh_kung;

% Plot Step Responses
systems_step_responses(systems)


