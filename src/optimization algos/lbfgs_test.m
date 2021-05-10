% %% FSOLVE LM ??
%% Params
mu = opt.g;
R_lm0 = Rini;
% LM_obj = @(R_lm) slra_mex_obj('func', obj, R_lm);
% LM_obj_reg = @(th) slra_mex_obj('func', obj, th2R(th)) ...
%                          + mu * norm(C(th), 'fro') ^ 2;
th_format.active    = 0;
th_format.PhiS_mat  = phi * (s0 + pext(tts + 1));
th_format.psi       = psi;    

%% Problem Formulation
if th_format.active
    LM_obj_reg = @(th) varproFuncAndGrad(obj, th2R(th), mu, 1, th_format);
    x0 = R2th(R_lm0, phi * (s0 + pext(tts + 1)), psi, opt.R0);
else                                           % Reg = 0
    LM_obj_reg = @(R) varproFuncAndGrad(obj, R, mu, 0, th_format);
    x0 = R_lm0;
end

% problem_lm.objective = LM_obj;
problem_lm.objective = LM_obj_reg;
problem_lm.x0 = x0;


%% fminlbfgs - bfgs
options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','bfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',0, 'MaxIter', 50);

tic
[x2,fval2,exitflag,output,grad] = fminlbfgs(LM_obj_reg,R_lm0,options);
toc
%% fminlbfgs - bfgs - GradConstr = false 
options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','bfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false, 'MaxIter', 30);
tic
[x2,fval2,exitflag,output,grad] = fminlbfgs(LM_obj_reg,R_lm0,options);
toc

%% fminlbfgs - lbfgs - GradConstr = false 
options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','lbfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false, 'MaxIter', 30);
tic
[x2,fval2,exitflag,output,grad] = fminlbfgs(LM_obj_reg,R_lm0,options);
toc

%%
sample_divisor1 = 1;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor1),:),y0(1:ceil(length(y0)/sample_divisor1),:), r2ss(x2, m_in, ell));

