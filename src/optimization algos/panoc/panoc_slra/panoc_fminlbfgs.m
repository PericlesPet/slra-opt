% PANOC algorithm using the lbgfs and proximal gradient methods
% include the Matlab folder with all its subfolders to your path
if isCloseAll == 1
    close all
    clc
end
%maxComplexity       = 15;
number_of_stepsBase = 300;
number_of_steps     = ceil(number_of_stepsBase * ...
    min((statsTable.complexities(cmplx_iter) / 500)^1.5,maxComplexity) / 100)*100;

printFreq       = 100;
max_number_of_steps_backtracking=8;
dimensions      = size(Rini);
dimension       = dimensions(1) * dimensions(2);
x_steps         = zeros(dimension,  number_of_steps);
grad_sizes      = zeros(3,          number_of_steps);
f_evals         = zeros(1,          number_of_steps);
f_lbfgs_evals   = zeros(1,          number_of_steps);
f_prox_evals    = zeros(1,          number_of_steps);
bcktrckingSteps = zeros(1,          number_of_steps);
taus            = zeros(1,number_of_steps);
gammas          = zeros(1,number_of_steps);
condition_array = zeros(max_number_of_steps_backtracking, number_of_steps, 5);
x0 = Rini(:);
% x0 = reshape(ss2r(sysh_kung), dimension, 1);
    % variables proximal gradient descent
beta_safety_value=0.05; % safety constant

g = @(x) g_indicator(reshape(x,dimensions));
proxg = @(x) prox_g_indicator(reshape(x,dimensions));

if ~exist('obj'), obj = slra_mex_obj('new', p, s, r); end    
mu = opt.g;
fcn_mode = 're';
[f, df] = SLRA_FuncAndGrad_handles(obj, mu, dimensions, fcn_mode);
LM_obj_reg = @(R) SLRA_FuncAndGrad_vals(obj, mu, dimensions, R, fcn_mode);

    % get starting value gamma
lipschitz_constant= estimate_lipschitz( df,x0 );
gamma0 = (1-beta_safety_value)/lipschitz_constant;
gamma=gamma0;
% gamma = 0.00001;
    % variables lbgfs
buffer_size=50; % buffer_size

    % internal buffers used
alpha=zeros(1,buffer_size);
beta=zeros(1,buffer_size);

    % ITERATE
    % About Display: 'off', 'iter', 'final', 'plot'
dsplLvl = 'off';
options = struct('GradObj','on','Display',dsplLvl,'LargeScale','off','HessUpdate','lbfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false, 'MaxIter', number_of_steps);
[exitflag, grad, data, optim] = fminlbfgs_data_init(LM_obj_reg,x0,options);
 
tau = 1;
bcktrckingSteps = zeros(1,          number_of_steps);
index           = 1;
clear panocFminlbfgsData
nonlcon = @(x) myConfungrad(x, constraints.ce_rri2, constraints.dce_rri2);

x=x0;   % set the starting point, and iterate
for interation_index=1:number_of_steps
    tic
    [x_prox_grad_descent,new_gamma] = prox_grad_descent_step( x,gamma,beta_safety_value,proxg,f,df );
    gamma=new_gamma;
    gammas(interation_index) = gamma;
    fval_prox  = f(x_prox_grad_descent);
    step_prox = (x-x_prox_grad_descent);
  
    if ~exist('x_lbfgs'), x_prev = x0(:); else, x_prev = data.xInitial; end
    [x_lbfgs,fval_lbfgs,exitflag,output, grad, data] = fminlbfgs_iteration_split1(LM_obj_reg, optim, data, exitflag, 0);

%%% CLASSIC LBFGS ITERATIONS
%     R= @(x) (1/gamma)*(x - proxg(x-df(x)*gamma));
%     [s_new,y_new,x_lbfgs] = lbfgs(interation_index,buffer_size,x,R,s,y);
%     s=s_new;y=y_new;
%%% Find the right convex combination trough backtracking
    tau=1;
    for i=1:max_number_of_steps_backtracking
        bcktrckingSteps(interation_index) = bcktrckingSteps(interation_index) + 1 ; 
        d=x_lbfgs(:)-x;
        potential_x=x + 1 * (-(1-tau)*step_prox + tau*d);
        sigma=beta_safety_value/(4*gamma);

    % The "Conditions" Are the Left and Right conditions for the backtracking step. 
        for get_conditions = 1  %This "for loop" is for Code Folding
            tau_condition = 2;
            %%% VARIOUS CONDITIONS FOR TAU
            switch tau_condition
                case 1      %%% NORMAL ONE      --- NOTE: CONDITION2 IS MAX(-0.1 * CONDITION1, ...)  
                    condition1           = FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg)  ;
                    condition2_first     = FBE( x,gamma,beta_safety_value,f,df,g,proxg);
                    condition2_second    = max(-0.1 * condition1, -sigma*norm(step_prox/gamma,2)^2);
                case 2      %%% BASED ON LBFGS AND PROX VALUES
%                     condition1           = FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg)  ;
                    condition1           = fval_lbfgs;
                    condition2_first     = fval_prox;
                    condition2_second    = 0.0001;       %max(-0.1 * condition1, -sigma*norm(step_prox/gamma,2)^2);                    
                case 3      %%% BASED ON SCALED FBE WITHOUT SECOND TERM (IF GAMMA TOO SMALL)
                    condition1           = FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg)  ;
                    condition2_first     = 1.1 * condition1;
                    condition2_second    = 0.0001; %max(-0.1 * condition1, -sigma*norm(step_prox/gamma,2)^2);
                otherwise
            end
            condition2  = condition2_first + condition2_second;            
            condition_t = toc;
            conditionz  = [condition1; ...
                         condition2_first; condition2_second; ... 
                         condition2; 1];
            condition_array(i , interation_index, :) = conditionz;
        end     

        if(condition1 <= condition2)            
            break; % if this is statified stop right away
        else
            tau=tau/2;
        end
    end
    
    
    
   
    x                               = potential_x;
%     x_steps(:,interation_index)     = x;
%     f_evals(interation_index)       = f(x);
%     f_lbfgs_evals(interation_index) = fval_lbfgs;
%     f_lbfgs_evals(interation_index) = f(x_lbfgs);
%     f_prox_evals(interation_index)  = fval_prox;
%     x_actuals(:, interation_index)  = x;
    taus(interation_index)          = tau;
    grad_sizes(:, interation_index) = [norm(x); norm(step_prox); norm(d)];

    updateProxFminLBFGS = 1;
    if updateProxFminLBFGS
            data.xInitial = potential_x;
            data.dir = (data.xInitial - x_prev)/data.alpha;
            data.alpha;
            % UPDATE HESSIAN VIA FminLBFGS    
            [x_lbfgs,fval_lbfgs,exitflag,output, grad, data] = ...
                fminlbfgs_iteration_split1(LM_obj_reg, optim, data, exitflag, 1);
    end            

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LOG DATA
    if exist('panocFminlbfgsData'), t_stamp = panocFminlbfgsData.t_stamps(end) + toc; ...
    else, t_stamp = toc; end
    R_panocfminlbfgs = reshape(x, size(slradata.Rini));     
    [~, M0] = compare(iddata(slradata.y0, slradata.u0), ...
        idss(r2ss(R_panocfminlbfgs, slradata.m_in, slradata.ell))); 

    panocFminlbfgsData.Xs(:, index)               = potential_x(:);
    panocFminlbfgsData.Xs_lbfgs(:, index)         = x_lbfgs(:);
    panocFminlbfgsData.Xs_proxgd(:, index)        = x_prox_grad_descent(:);

    panocFminlbfgsData.fvals(index)               = f(x); 
    panocFminlbfgsData.fvals_lbfgs(index)         = f(x_lbfgs); 
    panocFminlbfgsData.fvals_prox(index)          = f(x_prox_grad_descent); 
    panocFminlbfgsData.f_slra_val(index)          = slra_mex_obj('func', obj, R_panocfminlbfgs);

    [~, panocFminlbfgsData.CE(:,index), ~, ~]     = nonlcon(x(:));
    panocFminlbfgsData.M0(:,index)                = M0;
    panocFminlbfgsData.t_stamps(index)            = t_stamp;
    %     fminconData.M0(index)               = mean(M0);
    %         fminconData.CE(:,index)           = ce(x);    
    %         fminconData.CE(:,index)           = ce(x0); 

    index = index + 1;         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if mod(interation_index, printFreq) == 0
        fprintf('f(%d) = %f  -  t = %f  -  Accuracy = %f\n', ...
            interation_index, panocFminlbfgsData.fvals(interation_index),t_stamp, mean(M0)); end

end

if plotPANOC
    figure
    subplot(2,1,1)
    plot(panocFminlbfgsData.t_stamps,panocFminlbfgsData.fvals, 'k:', ...
        'Marker','*', 'MarkerSize', 3, 'MarkerIndices', 1:10:length(panocFminlbfgsData.t_stamps))
    hold on
    plot(panocFminlbfgsData.t_stamps,panocFminlbfgsData.fvals_lbfgs, 'r--')
    plot(panocFminlbfgsData.t_stamps,panocFminlbfgsData.fvals_prox, 'b--')
    legend('fvals', 'fvals_{lbfgs}', 'fvals_{prox}')
    ylim([0.9*f(R_true) max(panocFminlbfgsData.fvals)])
    title('F Evaluations')


    subplot(2,1,2)

    if ~isAccSemilog
        plot(panocFminlbfgsData.t_stamps,max(mean(panocFminlbfgsData.M0), 0))
    else
        semilogy(panocFminlbfgsData.t_stamps, ...
            max(mean(panocFminlbfgsData.M0(:,:)), ...
                1./abs(mean(panocFminlbfgsData.M0(:,:)))))

        ylim([0 100])
        % semilogy(panocFminlbfgsData.t_stamps, mean(panocFminlbfgsData.M0))
    end
    title('Mean Accuracy')
    suptitle('PANOC_{FMINLBFGS} Figures')
end
fprintf("Time Elapsed on PANOC_FMINLBFGS : %.3f\n", t_stamp);
R_panocfminlbfgs = reshape(x_steps(:,end), dimensions);

% [condition_xIndices, fbeLowCond, fbeHighCond] = FBE_ConditionProcessing(condition_array, 1, 4);
[condition_xIndices, fbeLowCond, fbeHighCond, fbeHigherCond] = FBE_ConditionProcessing(condition_array, 1, 4, 2);
% stiefConstraint(R_panoc)
% f_evals
% f_lbfgs_evals
% bcktrckingSteps

% %% plot the convergence rate
% figure(1);clf;
% subplot(2,1,1)
% plot(taus);
% 
% % figure(2);clf;
% subplot(2,1,2)
% plot(f_evals)    
% % % %%
% % % subplot_x = 2;
% % % subplot_y = 1;
% % % 
% % % figure;
% % % for subplot_1 = 1
% % %     subplot(subplot_x ,subplot_y, subplot_1)
% % %     plot(f_evals, 'ko', 'MarkerSize', 2)
% % %     hold on
% % %     plot(f_lbfgs_evals, '-.')
% % %     plot(f_prox_evals, '--')
% % %     plot(f_opt*ones(size(f_lbfgs_evals)))
% % %     legend('f evals', 'f lbfgs evals', 'f prox evals', 'optimal f (ground truth)')
% % % end
% % % 
% % % for subplot_2 = 2
% % %     subplot(subplot_x, subplot_y, subplot_2)
% % %     plot(f_evals, 'ko', 'MarkerSize', 2)
% % %     hold on
% % %     plot(f_lbfgs_evals, '-.')
% % %     plot(f_prox_evals, '--')
% % % %     plot(condition_xIndices, fbeLowCond, 'b', condition_xIndices, fbeHighCond, 'r')
% % %     plot(condition_xIndices, fbeLowCond, 'b', condition_xIndices, fbeHighCond, 'r', condition_xIndices, fbeHigherCond, '--k')
% % %     plot(f_opt*ones(size(f_lbfgs_evals)))
% % %     legend('f evals', 'f lbfgs evals', 'f prox evals', 'fbe low', 'fbe high', 'fbe higher', 'optimal f (ground truth)')
% % % end
% % % 
% % % %
% % % figure
% % % for subplot_3 = 1
% % %     subplot(subplot_x, subplot_y, subplot_3)
% % %     hold on
% % %     plot(gammas)
% % %     legend('gammas')
% % % end
% % % 
% % % for subplot_4 = 2
% % %     subplot(subplot_x, subplot_y, subplot_4)
% % %     plot(f_evals, 'ko', 'MarkerSize', 2)
% % %     hold on
% % % %     plot(condition_xIndices, fbeLowCond, 'b', condition_xIndices, fbeHighCond, 'r')
% % %     plot(condition_xIndices, fbeLowCond, 'b', condition_xIndices, fbeHighCond, 'r', condition_xIndices, fbeHigherCond, '--k')
% % %     legend('f evals', 'fbe low', 'fbe high', 'fbe higher')
% % % end
% % % 
% % % 
% % % 
% % % %%
% % % figure;
% % % plot(f_evals, 'o', 'MarkerSize', 2, 'MarkerFaceColor','b')
% % % hold on
% % % plot(condition_array(1,:))
% % % plot(condition_array(2,:))
% % % legend('f evals', 'fbe 1', 'fbe 2')
% % % % f(Rh)
% % % 
% % % 
% % % % figure(2);clf;
% plot(x_steps(1,2:end)./x_steps(1,1:end-1));