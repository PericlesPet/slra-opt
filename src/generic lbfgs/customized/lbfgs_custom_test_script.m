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

%%
fminlbfgs_iterations = 25;
x_fminlbfgs_steps = zeros(dimension, fminlbfgs_iterations);
f_fminlbfgs_vals  = zeros(1, fminlbfgs_iterations);
test_fs           = zeros(1, fminlbfgs_iterations);
%% fminlbfgs - lbfgs - GradConstr = false 
options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','lbfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false, 'MaxIter', 30);
tic

[exitflag, grad, data, optim] = fminlbfgs_data_init(LM_obj_reg,R_lm0,options);

% [x_next,fval_next,exitflag,output,grad, data] = fminlbfgs_iteration(LM_obj_reg,optim, data, exitflag);

noiseLevel = 0.001;
for i = 1:fminlbfgs_iterations
    
    [x_next,fval_next,exitflag,output,grad, data] = fminlbfgs_iteration_try(LM_obj_reg, optim, data, exitflag);
    x_fminlbfgs_steps(:, i) = x_next(:);
    f_fminlbfgs_vals(i) = fval_next;

%     data.xInitial = data.xInitial + noiseLevel * randn(size(data.xInitial));
    data.xInitial = x_next + noiseLevel * randn(size(data.xInitial));
    data.dir = (data.xInitial - x_fminlbfgs_steps(:,i))/data.alpha;
    [data,fval,grad]=gradient_function(data.xInitial,LM_obj_reg, data, optim);
    data.gradient=grad;
    
end
toc

for i = 1:fminlbfgs_iterations
    
    test_fs(i) = f(x_fminlbfgs_steps(:, i));
    
end
% %%
cheeeek = [f_fminlbfgs_vals; test_fs]
