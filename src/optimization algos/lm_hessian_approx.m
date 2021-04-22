clear logdata;
Rin = Rini;
lambda = 10^6;
for ii = 1:300
%%
f = slra_mex_obj('func', obj, Rin);
g = slra_mex_obj('grad', obj, Rin);
residual_R = Rin*Rin' ;
%%
hess_J_approx1 = g'*g;
sigma1 = -(hess_J_approx1 + lambda*eye(size(hess_J_approx1)))*g';

Rin = Rin + 0.01*sigma1' / norm(sigma1) ;
f_new = slra_mex_obj('func', obj, Rin);
%%
% hess_J_approx2 = g*g';
% sigma2 = -(hess_J_approx2 + lambda*eye(size(hess_J_approx2)))*g;
% 
% % f_new = slra_mex_obj('func', obj, Rin)
% Rin = Rin + 0.01*sigma2 / norm(sigma2) ;
% f_new = slra_mex_obj('func', obj, Rin);
%%
if (f_new <= f)
    lambda = lambda*0.5;
else
    lambda = lambda * 2;
end
%%
logdata(ii).f           = f_new;
% logdata(ii).f_reg       = f_reg;
logdata(ii).g           = g;
% logdata(ii).g_reg       = g_reg;
logdata(ii).residual_R  = residual_R;
logdata(ii).Rin         = Rin;
% logdata(ii).R_opt       = R - Rin;
% logdata(ii).dir_diff    = sign(R-Rin) - sign(sigma);
logdata(ii).lambda      = lambda;
%%
end
%
f_new
f_opt = min([logdata(:).f])


%% FSOLVE LM ??
%Params
mu = opt.g;
R_lm0 = Rini;
% LM_obj = @(R_lm) slra_mex_obj('func', obj, R_lm);
% LM_obj_reg = @(th) slra_mex_obj('func', obj, th2R(th)) ...
%                          + mu * norm(C(th), 'fro') ^ 2;
th_format.active    = 0;
th_format.PhiS_mat  = phi * (s0 + pext(tts + 1));
th_format.psi       = psi;    
        
%%
if th_format.active 
    
    LM_obj_reg = @(th) varproFuncAndGrad(obj, th2R(th), mu, 1, th_format);
    th_lm0 = R2th(R_lm0, phi * (s0 + pext(tts + 1)), psi, opt.R0);
    x0 = th_lm0;
else
    LM_obj_reg = @(R) varproFuncAndGrad(obj, R, mu, 1);
    x0 = R_lm0;
end

problem_lm.objective = LM_obj_reg;
problem_lm.x0 = x0;

% th_g = R2th(bb, phi * (s0 + pext(tts + 1)), psi, opt.R0);

%%
problem_lm.solver = 'fsolve';
problem_lm.options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter',... 
    'MaxFunctionEvaluations',10000,'MaxIterations',1500,'SpecifyObjectiveGradient',false);
problem_lm.options.TypicalX = ones(size(x0));

xopt = fsolve(problem_lm);

%%

%%
R_lmt = th2R(th_lmt);
f = slra_mex_obj('func', obj, R_lmt)
%%
f_opt = slra_mex_obj('func', obj, Ropt);
f_lm_init = slra_mex_obj('func', obj, R_lm0);
f_lm_opt = slra_mex_obj('func', obj, R_lmt);
[f_opt f_lm_init f_lm_opt]

residual_lmt = R_lmt*R_lmt';
residual_lm0 = R_lm0*R_lm0';
residual_opt = Ropt*Ropt';
[residual_opt residual_lm0 residual_lmt]


