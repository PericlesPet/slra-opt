%% simulation parameters
clear all, randn('seed', 0), rand('seed', 0), warning off
m = 2;          % Inputs
p = 2;          % Outputs
ell = 2;        % l time-horizon / dynamics                                                                                                  
n = ell * p;    % STATES -> outputs (X) * dynamics
q = m + p;      % w Vector size (INPUTS + OUTPUTS)

% %% TOTAL MATRIX ELEMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    n = l*p   STATES           %%%
%%%%    m         INPUTS           %%%
%%%%    p         OUTPUTS          %%%
%%%%    A -> n * n                 %%%
%%%%    B -> m * n                 %%%
%%%%    C -> n * p                 %%%
%%%%    D -> m * p                 %%%
%%%%    TOTAL: (n+m)*(p+n)         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total_elem = (n+m)*(p+n);
% %% SAMPLES, NOISE
T = 50*total_elem;     % Samples
s = 0.10;     % Noise Variation
opt_oe.exct = 1:m;      % fixed inputs = output error identification
opt_eiv.exct = [];      % errors-in-variables setup 
%% generate data

% Generate a discrete LTI system with:
% n STATES --- p OUTPUTS --- m INPUTS.
sys0 = drss(n, p, m);        % random "true" system

% RANDOM INPUT AND X_INIT
u0 = rand(T, m);        
xini0 = rand(n, 1);          % random "true" trajectory

% Y 
y0 = lsim(sys0, u0, [], xini0); 
yt = randn(T, p); 
y = y0 + s * norm(y0) * yt / norm(yt); % output error
u = u0;       

w0 = [u0 y0]; 
w = [u y];

%% compare ident
% profile on
% profiler_data = profile('info');
%   GET IDENT SOLUTION
tic, sysh_ident = ident(w, m, ell, opt_oe); t_ident = toc;

%   GET KUNG REALIZATION SOLUTION
tic, sysh_kung = w2h2ss(w, m, n); t_kung = toc;



%%
sys_comparison(u, y, sysh_ident, t_ident)
sys_comparison(u, y, sysh_kung, t_ident)

%% STEP RESPONSES
% Gather All Systems Into a Struct
systems.sys0 = sys0;
systems.sys1 = sysh_ident;
systems.sys2 = sysh_kung;

% Plot Step Responses
systems_step_responses(systems)


