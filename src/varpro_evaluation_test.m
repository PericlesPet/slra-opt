% Perform up to "Data Generation" of slra_sw_example
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
%%
% SLRA_MEX_OBJ 
% opt = rmfield(opt, 'solver');
obj = slra_mex_obj('new', p, s, r);

%%
% f1 = slra_mex_obj('func', 

tic, [ph, info] = slra_mex_obj('optimize', obj, opt); t_slra = toc

% [ph, info] = slra(w2p(w), s, r, opt); 
%%
[ph3, info3, wtferdata] = slra_ext(s2s(s, np), p, r, s.w, opt.Rini, s.phi, opt.psi, opt);

% %%
% prob.objective = @(th) Mslra_ext(th2R(th), tts, p, [], bfs, phi, s0);
% 
%                         Mslra_ext(R, tts, p, w, bfs, phi, s0)
%%
if exist('sysh_ident', 'var'), R_ident = ss2r(sysh_ident); end
if exist('sysh_kung', 'var'), R_kung = ss2r(sysh_kung); end
if exist('sys0', 'var'), R_opt = ss2r(sys0); end
Rh = info.Rh;
Rin = Rini;
%%
% Gradient Stepsize Parameter
gamma = 0.01;
Rin = Rini;
reg = exist('opt') && isfield(opt, 'method') && strcmp(opt.method, 'reg');
maxIter = 2000;
reg = 0;     % Use regularization
extCond = 1; % Use Exit condition to stop iterations
Ropt = ss2r(sys0);
tic
% [logdata, data_opt, f_log] = myGDesc(Ropt, maxIter, gamma, reg, opt, obj, R);
[logdata, data_opt, f_log, minf_log] = myGDesc(Rin, maxIter, gamma, reg, opt, obj, Rh, extCond);
% [logdata, data_opt, f_log, minf_log] = myProjGDesc(Rin, maxIter, gamma, reg, opt, obj, Rh, extCond);
gd_time = toc;

RmyOpt = data_opt.Rin;
data_opt

%%
% RmyOpt * RmyOpt'
% g = @(x) g_1(x,w);
f = @(R) slra_mex_obj('func', obj, R);

% f = slra_mex_obj('func', obj, R_try)
% R_try = stiefel_proj(Rini, RmyOpt);
% R_try * R_try'
%%
plot(f_log)
%%
numberz = 800;
(f_log(numberz +100)/f_log(numberz )-1)*100

%%
% sys_comparison(u0,y0, r2ss(R_try2, m_in, ell));
sample_divisor1 = 1;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor1),:),y0(1:ceil(length(y0)/sample_divisor1),:), r2ss(R_ident, m_in, ell));
sample_divisor2 = 1;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor2),:),y0(1:ceil(length(y0)/sample_divisor2),:), r2ss(RmyOpt, m_in, ell));
sample_divisor3 = 1;
sys_comparison(u(1:ceil(length(u)/sample_divisor3),:),y(1:ceil(length(y)/sample_divisor3),:), r2ss(Ropt, m_in, ell));
sample_divisor4 = 1;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor4),:),y0(1:ceil(length(y0)/sample_divisor4),:), r2ss(Rini, m_in, ell));
%%

%%
sample_divisor4 = 1;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor4),:),y0(1:ceil(length(y0)/sample_divisor4),:), r2ss(x, m_in, ell));
sys_comparison(u0(1:ceil(length(u0)/sample_divisor4),:),y0(1:ceil(length(y0)/sample_divisor4),:), r2ss(x2, m_in, ell));

%%
% So far:
% R          : slra_mex_obj('optimize') 
% RmyOpt     : optimal 
% Rini       : initial approximation
% sys0       : The Actual random system (ground-truth)
% sysh_ident : result of ident function
% sysh_kung  : result of Kung Realization

is_R = 1; 
if is_R  %Use SS format
    Rs.R = ss2r(sys0);
    Rs.Rh = Rh;    
    if ~exist('sysh_kung') Rs.Rkung = Rini; else Rs.Rkung = ss2r(sysh_kung); end
    Rs.RmyOpt = RmyOpt;    

    SSs = [];
else % Use SS format
    SSs.sys0        = sys0;
    if exist('sysh_ident') SSs.sysh_ident  = sysh_ident; else SSs.sysh_ident = r2ss(Rh); end
    if exist('sysh_kung') SSs.sysh_kung = sysh_kung; else SSs.sysh_kung = r2ss(Rini); end
    SSs.sysh_myOpt  = r2ss(RmyOpt);

    RRs = [];
end
%%

M_S = phi * (s0 + pext(tts + 1));
R_stats = R_comparison(Rs,SSs, is_R, C, obj, m_in, ell, u0, y0, th2R, M_S, psi);


%% DELETE OBJ
slra_mex_obj('delete', obj);
