clear logdata;
Rin = Rini;
lambda = 10^4;
for ii = 1:10
%%
f = slra_mex_obj('func', obj, Rin);
g = slra_mex_obj('grad', obj, Rin);
residual_R = Rin*Rin' ;
%%
hess_J_approx = g'*g;
sigma = -(hess_J_approx + lambda*eye(size(hess_J_approx)))*g';

Rin = Rin + 0.001*sigma' ;
f_new = slra_mex_obj('func', obj, Rin);
if (f_new <= f)
    lambda = lambda*0.5;
else
    lambda = lambda * 2;
end
%%
logdata(ii).f           = f;
% logdata(ii).f_reg       = f_reg;
logdata(ii).g           = g;
% logdata(ii).g_reg       = g_reg;
logdata(ii).residual_R  = residual_R;
logdata(ii).Rin         = Rin;
logdata(ii).R_opt       = R - Rin;
logdata(ii).dir_diff    = sign(R-Rin) - sign(sigma);
%%
end




%% FSOLVE LM ??

R_lm0 = Rin;
% R_lm0 = ones(size(Rin));
LM_obj = @(R_lm) slra_mex_obj('func', obj, R_lm);
LM_obj_reg = @(th) slra_mex_obj('func', obj, th2R(th)) ...
                         + mu * norm(C(th), 'fro') ^ 2;

% problem_lm.objective = LM_obj;
problem_lm.objective = LM_obj_reg;
% problem_lm.x0 = R_lm0;
problem_lm.x0 = R2th(R_lm0, phi * (s0 + pext(tts + 1)), psi, opt.R0); 

problem_lm.solver = 'fsolve';
problem_lm.options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','iter',... 
    'MaxFunctionEvaluations',4500,'MaxIterations',1500);

%%
R_lmt = fsolve(problem_lm)

%%
f_opt = slra_mex_obj('func', obj, R);
f_lm_init = slra_mex_obj('func', obj, R_lm0);
f_lm_opt = slra_mex_obj('func', obj, R_lmt);
[f_opt f_lm_init f_lm_opt]

residual_lmt = R_lmt*R_lmt';
residual_lm0 = R_lm0*R_lm0';
residual_opt = R*R';
[residual_opt residual_lm0 residual_lmt]


