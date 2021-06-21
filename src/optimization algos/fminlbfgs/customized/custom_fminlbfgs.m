%% Objective Params
close all
clc
mu = opt.g;
R_lm0 = Rini;

th_format.active    = 0;
th_format.PhiS_mat  = phi * (s0 + pext(tts + 1));
th_format.psi       = psi;    

    % Problem Formulation
if th_format.active
    LM_obj_reg = @(th) varproFuncAndGrad(obj, th2R(th), mu, 1, th_format);
    x0 = R2th(R_lm0, phi * (s0 + pext(tts + 1)), psi, opt.R0);
else                                           % Reg = 0
    LM_obj_reg = @(R) SLRA_FuncAndGrad_vals(obj, mu, dimensions, R, 'reg');
    x0 = R_lm0;
end

problem_lm.objective = LM_obj_reg;
problem_lm.x0 = x0;

% %% FminLBFGS Params
fminlbfgs_iterations    = 400;
pcntImproveThresh       = 1.0;
selectUpdate            = 1;
noiseLevel              = 0.00000;
displayComparisons      = 0;

x_fminlbfgs_steps       = zeros(dimension, fminlbfgs_iterations);
f_fminlbfgs_vals        = zeros(1, fminlbfgs_iterations);
f_fminlbfgs_prox_vals   = zeros(1, fminlbfgs_iterations);
test_fs                 = zeros(1, fminlbfgs_iterations);


    % ITERATE
    % fminlbfgs - lbfgs - GradConstr = false 
clear fminlbfgsData    
    % About Display: 'off', 'iter', 'final', 'plot'
dsplLvl = 'off';
options = struct('GradObj','on','Display',dsplLvl,'LargeScale','off','HessUpdate','lbfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false, 'MaxIter', fminlbfgs_iterations);

% [exitflag, grad, data, optim] = fminlbfgs_data_init(LM_obj_reg,x0,options);
[exitflag, grad, data, optim] = fminlbfgs_data_init(LM_obj_reg,R_lm0,options);
% [exitflag, grad, data, optim] = fminlbfgs_data_init(LM_obj_reg,Rkung,options);
index = 1;
for i = 1:fminlbfgs_iterations

    tic
    % STORE PREVIOUS X VALUE, X_PREV
    if ~exist('x_lbfgs'), x_prev = R_lm0(:); else, x_prev = data.xInitial; end
    if ~isfield(data, 'dir'), dir_prev = zeros(size(x_prev)); else, dir_prev = data.dir; end
    
    % UPDATE X VIA FminLBFGS
    [x_lbfgs,fval_lbfgs,exitflag,output, grad, data] = ...
        fminlbfgs_iteration_split1(LM_obj_reg, optim, data, exitflag, 0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fminlbfgsData.Xs_lbfgs(:, index)        = x_lbfgs(:);
    fminlbfgsData.fvals_lbfgs(index)    = fval_lbfgs; 

%     x_fminlbfgs_steps(:, i) = x_lbfgs(:);
%     f_fminlbfgs_vals(i) = fval_lbfgs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Update X, Direction, 
    xInitial_prev = data.xInitial;
    if selectUpdate == 0
        fprintf('Update fminlbfgs with random noise (to check robustness)\n');
        data.xInitial = x_lbfgs(:) + noiseLevel * randn(size(data.xInitial));
        data.dir = (data.xInitial - x_prev);        
        step_prox = x_lbfgs(:) - data.xInitial(:);
        potential_x = data.xInitial;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fminlbfgsData.fvals_prox(index) = f(potential_x);
        fminlbfgsData.Xs_prox(:, index)          = potential_x(:);
        %         f_fminlbfgs_prox_vals(i) = f(potential_x);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    elseif selectUpdate == 1 % Update via PROX
        %         fprintf('Update fminlbfgs with proximal step\n');
        [x_prox_grad_descent,new_gamma] = prox_grad_descent_step( x_lbfgs(:),gamma,beta_safety_value,proxg,f,df );
        gamma=new_gamma;
        step_prox = (x_lbfgs(:)-x_prox_grad_descent);
        potential_x = x_lbfgs(:) - step_prox;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fminlbfgsData.fvals_prox(index) = f(potential_x);
        fminlbfgsData.Xs_prox(:, index)          = potential_x(:);
%         f_fminlbfgs_prox_vals(i) = f(potential_x);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        prcnt_change = (fminlbfgsData.fvals_lbfgs(i) - fminlbfgsData.fvals_prox(i)) ./ (fminlbfgsData.fvals_lbfgs(i)) * 100 ;

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
    
    
    % Display some things for diagnostics
    if displayComparisons
        x_compare = [x_prev(1:5)'; xInitial_prev(1:5)' ; ...
                    data.xInitial(1:5)' ; potential_x(1:5)']
        dir_next = (data.xInitial - fminlbfgsData.Xs_lbfgs(:,i));
        dir_compare = [dir_prev(1:5)' ; ...
                    (data.xInitial(1:5) - x_prev(1:5))' ; ...
                    (1/data.alpha)*(data.xInitial(1:5) - x_prev(1:5))' ; ...
                    (dir_prev(1:5) - step_prox(1:5))' ; ...
                    dir_next(1:5)' ] ;
        grad_compare = [grad_prev(1:5)' ; data.gradient(1:5)'];
    end

    % UPDATE HESSIAN VIA FminLBFGS    
    [x_lbfgs,fval_lbfgs,exitflag,output, grad, data] = ...
        fminlbfgs_iteration_split1(LM_obj_reg, optim, data, exitflag, 1);

    % LOG DATA
    if isfield(fminlbfgsData,'t_stamps'), t_stamp = fminlbfgsData.t_stamps(end) + toc; ...
    else, t_stamp = toc; end
    
    R_fminlbfgs = x_lbfgs;     
    [~, M0] = compare(iddata(slradata.y0, slradata.u0), ...
        idss(r2ss(R_fminlbfgs, slradata.m_in, slradata.ell))); 
    fminlbfgsData.Xs(:, index)               = x_lbfgs(:);
    fminlbfgsData.fvals(index)               = fval_lbfgs; 
    if selectUpdate ~= 1
        fminlbfgsData.fvals_prox(index) = 0;
    end
    %     fminconData.M0(index)               = mean(M0);
    fminlbfgsData.M0(:,index)                = M0;
    fminlbfgsData.f_slra_val(index)          = slra_mex_obj('func', obj, R_fminlbfgs);
    %     fminconData.CE(:,index)           = ce(x);    
    %     fminconData.CE(:,index)           = ce(x0); 
    [~, fminlbfgsData.CE(:,index), ~, ~]     = nonlcon(x_lbfgs(:));
    fminlbfgsData.t_stamps(index)         = t_stamp;

    index = index + 1;         

    
    
    if mod(i, 20) == 0
        fprintf('f(%d) = %f  -  t = %f  -  Accuracy = %f\n', ...
            i, panocLbfgsData.fvals(i),t_stamp, mean(M0)); end


    
end



check1 = [fminlbfgsData.fvals_lbfgs; ...
        fminlbfgsData.fvals_prox];

check2 = [fminlbfgsData.fvals_lbfgs; ...
        fminlbfgsData.fvals_prox; ... 
        ((fminlbfgsData.fvals_lbfgs(:) - fminlbfgsData.fvals_prox(:)) ...
        ./ (fminlbfgsData.fvals_lbfgs(:)) * 100)'];

figure
subplot(2,1,1)
plot(fminlbfgsData.t_stamps,fminlbfgsData.fvals, 'k:', ...
    'Marker','*', 'MarkerSize', 3, 'MarkerIndices', 1:10:length(fminlbfgsData.t_stamps))
hold on
plot(fminlbfgsData.t_stamps,fminlbfgsData.fvals_lbfgs, 'r--')
plot(fminlbfgsData.t_stamps,fminlbfgsData.fvals_prox, 'b--')
legend('fvals', 'fvals_{lbfgs}', 'fvals_{prox}')
ylim([0.9*f(R_true) max(fminlbfgsData.fvals)])
title('F Evaluations')


subplot(2,1,2)
plot(fminlbfgsData.t_stamps,mean(fminlbfgsData.M0))
title('Mean Accuracy')

fprintf("Time Elapsed on PANOC_FMINLBFGS : %.3f\n", t_stamp);
R_panocfminlbfgs = reshape(x_steps(:,end), dimensions);
