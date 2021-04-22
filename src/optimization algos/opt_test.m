
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
LM_obj_reg = @(th) varproFuncAndGrad(obj, th2R(th), mu, 1, th_format);

% problem_lm.objective = LM_obj;
problem_lm.objective = LM_obj_reg;
th_lm0 = R2th(R_lm0, phi * (s0 + pext(tts + 1)), psi, opt.R0);

% th_g = R2th(bb, phi * (s0 + pext(tts + 1)), psi, opt.R0);

problem_lm.x0 = th_lm0;
%%



%%
options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','bfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',0, 'MaxIter', 40);
% LM_obj_reg = @(th) varproFuncAndGrad(obj, th2R(th), mu, 1, th_format);
% th_lm0 = R2th(R_lm0, phi * (s0 + pext(tts + 1)), psi, opt.R0);
LM_obj_reg = @(R) varproFuncAndGrad(obj, R, mu, 1, th_format);
R_lm0 = Rini;
% R_lm0 = Rin;
%%
tic
[x2,fval2,exitflag,output,grad] = fminlbfgs(LM_obj_reg,R_lm0,options);
toc
%%
options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','bfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false, 'MaxIter', 30);
tic
[x2,fval2,exitflag,output,grad] = fminlbfgs(LM_obj_reg,R_lm0,options);
toc
%%
tic
[x,fval] = fminunc(LM_obj_reg,R_lm0,options);
toc


%%
sample_divisor1 = 1;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor1),:),y0(1:ceil(length(y0)/sample_divisor1),:), r2ss(x2, m_in, ell));

