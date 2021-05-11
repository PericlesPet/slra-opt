%% Objective Params
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

%% FminLBFGS Params
fminlbfgs_iterations    = 50;
pcntImproveThresh       = 0.5;
useNoise                = 0;
noiseLevel              = 0.00001;
displayComparisons      = 0;

x_fminlbfgs_steps       = zeros(dimension, fminlbfgs_iterations);
f_fminlbfgs_vals        = zeros(1, fminlbfgs_iterations);
f_fminlbfgs_prox_vals   = zeros(1, fminlbfgs_iterations);
test_fs                 = zeros(1, fminlbfgs_iterations);

%% fminlbfgs - lbfgs - GradConstr = false 

options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','lbfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false, 'MaxIter', fminlbfgs_iterations);
tic

[exitflag, grad, data, optim] = fminlbfgs_data_init(LM_obj_reg,R_lm0,options);

for i = 1:fminlbfgs_iterations-2

%     STORE PREVIOUS X VALUE, X_PREV
    if ~exist('x_lbfgs'), x_prev = R_lm0(:); else, x_prev = data.xInitial; end
    if ~isfield(data, 'dir'), dir_prev = zeros(size(x_prev)); else, dir_prev = data.dir; end
    
% UPDATE X VIA FminLBFGS
    [x_lbfgs,fval_lbfgs,exitflag,output, grad, data] = fminlbfgs_iteration_split1(LM_obj_reg, optim, data, exitflag, 0);

    x_fminlbfgs_steps(:, i) = x_lbfgs(:);
    f_fminlbfgs_vals(i) = fval_lbfgs;
    

    if fval_lbfgs <= 72
        break;
    end 

% Update X, Direction, 
    xInitial_prev = data.xInitial;
    if useNoise
        data.xInitial = x_lbfgs(:) + noiseLevel * randn(size(data.xInitial));
        data.dir = (data.xInitial - x_prev);        
        step_prox = x_lbfgs(:) - data.xInitial(:);
    else % Update via PROX
        [x_prox_grad_descent,new_gamma] = prox_grad_descent_step( x_lbfgs(:),gamma,beta_safety_value,proxg,f,df );
        gamma=new_gamma;
        step_prox = (x_lbfgs(:)-x_prox_grad_descent);
        potential_x = x_lbfgs(:) - step_prox;
        f_fminlbfgs_prox_vals(i) = f(potential_x);
        prcnt_change = (f_fminlbfgs_vals(i) - f_fminlbfgs_prox_vals(i)) ./ (f_fminlbfgs_vals(i)) * 100 ;

        if prcnt_change >= pcntImproveThresh
            data.xInitial = potential_x;
            data.dir = (data.xInitial - x_prev)/data.alpha;
            data.alpha;
        end            
    end
    
    
% Update Gradient
    [data,fval,grad]=gradient_function(data.xInitial,LM_obj_reg, data, optim);
    grad_prev = data.gradient;
    data.gradient=grad;
    
    

    if displayComparisons
        x_compare = [x_prev(1:5)'; xInitial_prev(1:5)' ; data.xInitial(1:5)' ; potential_x(1:5)']
        dir_next = (data.xInitial - x_fminlbfgs_steps(:,i));
        dir_compare = [dir_prev(1:5)' ; ...
                    (data.xInitial(1:5) - x_prev(1:5))' ; ...
                    (1/data.alpha)*(data.xInitial(1:5) - x_prev(1:5))' ; ...
                    (dir_prev(1:5) - step_prox(1:5))' ; ...
                    dir_next(1:5)' ]
        grad_compare = [grad_prev(1:5)' ; data.gradient(1:5)']
    end

% UPDATE HESSIAN VIA FminLBFGS    
    [x_lbfgs,fval_lbfgs,exitflag,output, grad, data] = fminlbfgs_iteration_split1(LM_obj_reg, optim, data, exitflag, 1);

end
toc

check1 = [f_fminlbfgs_vals; ...
        f_fminlbfgs_prox_vals];

check2 = [f_fminlbfgs_vals; ...
        f_fminlbfgs_prox_vals; ... 
        ((f_fminlbfgs_vals(:) - f_fminlbfgs_prox_vals(:)) ./ (f_fminlbfgs_vals(:)) * 100)']

x_finalz = x_fminlbfgs_steps(:,i);
sys_comparison(u0,y0, r2ss(reshape(x_finalz, dimensions), m_in, ell));
