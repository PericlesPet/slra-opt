% Perform up to "Data Generation" of slra_sw_example

[w, s, r, opt, q, N, T] = sys2slra(w, m_in, ell, opt_oe);

p = w2p(w);
R = ss2r(sys0); 
np = length(p); 

[tts, p, r, s, w, Rini, phi ,psi, opt, th2R, C, s0, prob] = slra2slra_ext(p, s, r, opt);
  
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
gamma = 0.001;

reg = exist('opt') && isfield(opt, 'method') && strcmp(opt.method, 'reg');
maxIter = 500;
[logdata, data_opt] = myGDesc(Rin, maxIter, gamma, reg, opt, obj, R);
RmyOpt = data_opt.Rin;

data_opt
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
    if ~exist('sysh_kung') Rs.Rkung = R; else Rs.Rkung = ss2r(sysh_kung); end
    Rs.RmyOpt = RmyOpt;    

    SSs = [];
else % Use SS format
    SSs.sys0        = sys0;
    if exist('sysh_ident') SSs.sysh_ident  = sysh_ident; else SSs.sysh_ident = r2ss(R); end
    if exist('sysh_kung') SSs.sysh_kung = sysh_kung; else SSs.sysh_kung = r2ss(R); end
    SSs.sysh_myOpt  = r2ss(RmyOpt);

    RRs = [];
end
%%
R_stats = R_comparison(Rs,SSs, is_R, C, obj, m_in, ell, u0, y0);


%% DELETE OBJ
slra_mex_obj('delete', obj);
