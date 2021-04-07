%% simulation parameters
clear all, randn('seed', 0), rand('seed', 0), warning off
m = 3;          % Inputs
p = 3;          % Outputs
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
s = 0.001;     % Noise Variation
opt_oe.exct = 1:m;      % fixed inputs = output error identification
opt_eiv.exct = [];      % errors-in-variables setup 
%% generate data

% Generate a discrete LTI system with:
% n STATES --- p OUTPUTS --- m INPUTS.
sys0 = drss(n, p, m);        % random "true" system
u0 = rand(T, m); 
xini0 = rand(n, 1);          % random "true" trajectory

% Y 
y0 = lsim(sys0, u0, [], xini0); 
u = u0;       
yt = randn(T, p); 
y = y0 + s * norm(y0) * yt / norm(yt); % output error

w0 = [u0 y0]; 
w = [u y];

%% compare ident
profile on
tic, sysh_ident = ident(w, m, ell, opt_oe); t_ident = toc;
profiler_data = profile('info');
% [Yh, M] = compare(iddata(y, u), idss(sysh_ident), sysh_pem); 
figure
compare(iddata(y, u), idss(sysh_ident)); 
[Yh, M] = compare(iddata(y, u), idss(sysh_ident)); 
t_ident 
M

%% STEP RESPONSES
figure
t = 1:1:100;
step(sys0,t)
hold on
step(sysh_ident,t)
% step(sysh_pem)
% legend('real system', 'IDed system', 'PEM')
legend('real system', 'IDed system')
% step(sysh_pem)

