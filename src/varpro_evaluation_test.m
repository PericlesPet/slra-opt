% Perform up to "Data Generation" of slra_sw_example
tic
[w, s, r, opt, q, N, T] = sys2slra(w, m_in, ell, opt_mo);
toc
p = w2p(w);
R = ss2r(sys0); 
np = length(p); 
tic
[tts, p, r, s, w, Rini, phi ,psi, opt, th2R, C, s0, prob, pext] = slra2slra_ext(p, s, r, opt);
toc
%% SLRA_MEX_OBJ 

obj = slra_mex_obj('new', p, s, r);
%%
[ph, info] = slra_mex_obj('optimize', obj, opt);

if exist('sysh_ident', 'var'), R_ident = ss2r(sysh_ident); end
R = info.Rh;
% Rin = ones(size(R));
Rin = Rini;

%%
% Gradient Stepsize Parameter
gamma = 0.01;
Rin = Rini;
reg = exist('opt') && isfield(opt, 'method') && strcmp(opt.method, 'reg');
maxIter = 10000;
reg = 0;
Ropt = ss2r(sys0);
tic
[logdata, data_opt, f_log] = myGDesc(Ropt, maxIter, gamma, reg, opt, obj, R);
gd_time = toc;
% [logdata, data_opt] = myGDesc(Rin, maxIter, gamma, reg, opt, obj, R);

RmyOpt = data_opt.Rin;
data_opt

%% stiefel
R_try = stiefel_proj(Rini, RmyOpt);
%% stiefel
R_try2 = stiefel_proj(Rini, R_try);
%% stiefel
R_try2 = stiefel_proj(Rini, R_try2);
%% 
RmyOpt * RmyOpt'
R_try * R_try'
f = slra_mex_obj('func', obj, R_try)
%%
% sys_comparison(u0,y0, r2ss(R_try2, m_in, ell));
sample_divisor1 = 60;
sample_divisor2 = 30;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor1),:),y0(1:ceil(length(y0)/sample_divisor1),:), r2ss(R_try, m_in, ell));
sys_comparison(u0(1:ceil(length(u0)/sample_divisor2),:),y0(1:ceil(length(y0)/sample_divisor2),:), r2ss(RmyOpt, m_in, ell));
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
    Rs.Rh = R;    
    if ~exist('sysh_kung') Rs.Rkung = Rini; else Rs.Rkung = ss2r(sysh_kung); end
    Rs.RmyOpt = RmyOpt;    

    SSs = [];
else % Use SS format
    SSs.sys0        = sys0;
    if exist('sysh_ident') SSs.sysh_ident  = sysh_ident; else SSs.sysh_ident = r2ss(R); end
    if exist('sysh_kung') SSs.sysh_kung = sysh_kung; else SSs.sysh_kung = r2ss(Rini); end
    SSs.sysh_myOpt  = r2ss(RmyOpt);

    RRs = [];
end
%%

M_S = phi * (s0 + pext(tts + 1));
R_stats = R_comparison(Rs,SSs, is_R, C, obj, m_in, ell, u0, y0, th2R, M_S, psi);


%% DELETE OBJ
slra_mex_obj('delete', obj);
