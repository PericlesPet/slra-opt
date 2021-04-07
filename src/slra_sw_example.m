%% simulation parameters
clear all, randn('seed', 0), rand('seed', 0), warning off
m = 3;       % Inputs
p = 3;       % Outputs
ell = 3;     % l time-horizon / dynamics                                                                                                  
T = 300;     % Samples
s = 0.1;     % Noise Variation
opt_oe.exct = 1:m;      % fixed inputs = output error identification
opt_eiv.exct = [];      % errors-in-variables setup 

%% generate data
n = ell * p;    % STATES -> outputs (X) * dynamics
q = m + p;      
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

%% compare ident (and pem)
profile on


tic, sysh_ident = ident(w, m, ell, opt_oe); t_ident = toc;

profiler_data = profile('info')
t_ident
% tic, sysh_pem = pem(iddata(y, u), n); t_pem = toc;
% [Yh, M] = compare(iddata(y, u), idss(sysh_ident), sysh_pem); 
% % % figure
% % % compare(iddata(y, u), sysh_ident); 
% % % % ans = [t_ident; t_pem]
% % % ans = [t_ident]


%%
figure
step(sys0)
hold on
step(sysh_ident)
% step(sysh_pem)
% legend('real system', 'IDed system', 'PEM')
legend('real system', 'IDed system')
% step(sysh_pem)






%% simulate step response data
ys0 = step(sys0);       %Step Response of sys0
Ts = length(ys0);        
us = ones(Ts, m);       % Step Input
yt = randn(Ts, p); 
ys = ys0 + s * norm(ys0) * yt / norm(yt);       % Noisy

opt_oe.w0 = 0; opt_oe.sys0 = sys0; 
[sysh, info, wh] = ident([us ys], m, ell, opt_oe);

% PLOT
figure;
plot(1:Ts, ys0, '-r', 1:Ts, ys, ':k', 1:Ts, wh(:, end), '--b')
legend('true', 'noisy', 'approx.   .')
