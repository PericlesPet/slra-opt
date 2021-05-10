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

noiseLevel = 0.00001;
for i = 1:fminlbfgs_iterations

    if ~exist('x_next'), x_prev = R_lm0(:); else, x_prev = data.xInitial; end;
    [x_next,fval_next,exitflag,output, grad, data] = fminlbfgs_iteration_try(LM_obj_reg, optim, data, exitflag);
    x_fminlbfgs_steps(:, i) = x_next(:);
    f_fminlbfgs_vals(i) = fval_next;
    

%     data.xInitial = data.xInitial + noiseLevel * randn(size(data.xInitial));

% Update X via PROX
    [x_prox_grad_descent,new_gamma] = prox_grad_descent_step( x_next(:),gamma,beta_safety_value,proxg,f,df );
    gamma=new_gamma;
    step_prox = (x_next(:)-x_prox_grad_descent);
    potential_x = x_next(:) - step_prox;

% Update X
    xInitial_prev = data.xInitial;
    data.xInitial = x_next(:) + noiseLevel * randn(size(data.xInitial));
    x_compare = [x_prev'; xInitial_prev' ; data.xInitial' ; potential_x']
% Update Direction
    %     data.dir = (data.xInitial - x_fminlbfgs_steps(:,i))/data.alpha;
    dir_prev = data.dir;
    dir_next = (data.xInitial - x_fminlbfgs_steps(:,i));
    data.dir = (data.xInitial - x_prev);
    
    dir_compare = [dir_prev' ; data.dir' ; dir_next' ; step_prox']
    alignability = (dir_prev' * data.dir) / (norm(dir_prev)^2 )
    
% Update Gradient
    grad_prev = data.gradient;
    [data,fval,grad]=gradient_function(data.xInitial,LM_obj_reg, data, optim);
    data.gradient=grad;
    grad_compare = [grad_prev' ; data.gradient']

end
toc

for i = 1:fminlbfgs_iterations
    
    test_fs(i) = f(x_fminlbfgs_steps(:, i));
    
end
% %%
cheeeek = [f_fminlbfgs_vals; test_fs]



function [data,fval,grad]=gradient_function(x,funfcn, data, optim)
    % Call the error function for error (and gradient)
    if ( nargout <3 )
        timem=tic;   
        fval=funfcn(reshape(x,data.xsizes)); 
        data.timeExtern=data.timeExtern+toc(timem);
        data.funcCount=data.funcCount+1;
    else
        if(strcmp(optim.GradObj,'on'))
            timem=tic;    
            [fval, grad]=feval(funfcn,reshape(x,data.xsizes)); 
            data.timeExtern=data.timeExtern+toc(timem);
            data.funcCount=data.funcCount+1;
            data.gradCount=data.gradCount+1;
        else
            % Calculate gradient with forward difference if not provided by the function
            grad=zeros(length(x),1);
            fval=funfcn(reshape(x,data.xsizes));
            gstep=data.initialStepLength/1e6; 
            if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
            if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
            for i=1:length(x),
                x_temp=x; x_temp(i)=x_temp(i)+gstep;
                timem=tic;    
                [fval_g]=feval(funfcn,reshape(x_temp,data.xsizes)); data.funcCount=data.funcCount+1;
                data.timeExtern=data.timeExtern+toc(timem);
                grad(i)=(fval_g-fval)/gstep;
            end
        end
        grad=grad(:);
    end
end