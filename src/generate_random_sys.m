function [n, q, T, sys0, u0, y0, w0, u, y, w] = generate_random_sys(m, p, ell, s)
% INPUT ARGUMENTS
% m     : Inputs
% p     : Outputs
% ell   : l time-horizon / dynamics                                                                                                  
% s     : Noise Variation

% OUTPUT ARGUMENTS
% n,q           : States, Vector Size (Inputs + Outputs)
% T             : Samples
% sys0          : State Space Model of random system
% u0, y0, w0    : (Initial System)    u0 -> Input Vector , y0 -> Output Vector, w0 -> [u0 y0]
% u, y, w       : (With Added Noise)  u  -> Input Vector , y  -> Output Vector, w  -> [u y]


n = ell * p;    % STATES -> outputs (X) * dynamics
q = m + p;      % w Vector size (INPUTS + OUTPUTS)

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
total_elem = (n+m)*(p+n);
% %% SAMPLES, NOISE
T = 50*total_elem;     % Samples, based on system complexity
% T = 10;     % Samples

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


