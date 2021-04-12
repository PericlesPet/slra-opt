%% simulation parameters
% addpath(genpath('helper_functions'))
addpath(genpath('..\..\..\Matlab'))
clear all, randn('seed', 0), rand('seed', 0), warning off
m_in = 1;            % Inputs
p_out = 1;           % Outputs
ell = 2;             % l time-horizon / dynamics                                                                                                  
s_noise = 0.10;      % Noise Variation

%% GENERATE RANDOM SYSTEM
[n, q, T, sys0, u0, y0, w0, u, y, w] = generate_random_sys(m_in, p_out, ell, s_noise);

% SLRA Solver Options
opt_oe.exct = 1:m_in;      % fixed inputs = output error identification
opt_oe.solver = 'm';
% opt_oe.method = 'reg';
opt_eiv.exct = [];      % errors-in-variables setup 

%%
for not_happening = 1:-1:2 
not_happening    
n = ell * p_out;    % STATES -> outputs (X) * dynamics
q = m_in + p_out;      % w Vector size (INPUTS + OUTPUTS)

% %% TOTAL MATRIX ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    n = l*p      STATES        %%%
%%%%    m            INPUTS        %%%
%%%%    p            OUTPUTS       %%%
%%%%    A -> n * n   ELEMENTS      %%%
%%%%    B -> m * n   ELEMENTS      %%%
%%%%    C -> n * p   ELEMENTS      %%%
%%%%    D -> m * p   ELEMENTS      %%%
%%%%_______________________________%%%
%%%%  = (n+m)*(p+n) TOTAL ELEMENTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_elem = (n+m_in)*(p_out+n);
% %% SAMPLES, NOISE
T = 50*total_elem;     % Samples, based on system complexity
% T = 10;     % Samples
%% generate data

% Generate a discrete LTI system with:
% n STATES --- p OUTPUTS --- m INPUTS.
sys0 = drss(n, p_out, m_in);        % random "true" system

% RANDOM INPUT AND X_INIT
u0 = rand(T, m_in);        
xini0 = rand(n, 1);          % random "true" trajectory

% Y 
y0 = lsim(sys0, u0, [], xini0); 
yt = randn(T, p_out); 
y = y0 + s_noise * norm(y0) * yt / norm(yt); % output error
u = u0;       

w0 = [u0 y0]; 
w = [u y];

% w,m, ell, opt_oe, n, u, y, sys0

end
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


