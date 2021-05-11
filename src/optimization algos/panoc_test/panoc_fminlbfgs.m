% PANOC algorithm using the lbgfs and proximal gradient methods
% include the Matlab folder with all its subfolders to your path
number_of_steps = 20;
dimensions      = size(Rini);
dimension       = dimensions(1) * dimensions(2);
x_steps         = zeros(dimension,  number_of_steps);
grad_sizes      = zeros(3,          number_of_steps);
f_evals         = zeros(1,          number_of_steps);
f_lbfgs_evals   = zeros(1,          number_of_steps);
f_prox_evals    = zeros(1,          number_of_steps);
bcktrckingSteps = zeros(1,          number_of_steps);
taus            = zeros(1,number_of_steps);
condition_array = zeros(6,number_of_steps);
x0 = Rini(:);
% x0=[0.5;0.5];
% x0=100;
% %% variables proximal gradient descent
beta_safety_value=0.05; % safety constant

% dimensions = size(Rini); 
g = @(x) g_indicator(reshape(x,dimensions));
proxg = @(x) prox_g_indicator(reshape(x,dimensions));
% g = @(x) g_2(x);
% proxg = @(x) prox_g_2( x );

% degree_polynomial=20;
% f = @(x) sum(x.^degree_polynomial);
% df = @(x) (degree_polynomial)*x.^(degree_polynomial-1);

if ~exist('obj'), obj = slra_mex_obj('new', p, s, r); end    
% f   = @(x) slra_mex_obj('func', obj, reshape(x, dimensions));
% df  = @(x) slra_mex_obj('grad', obj, reshape(x, dimensions));
f   = @(x) slra_mex_obj('func', obj, reshape(x, dimensions));
df  = @(x) reshape(slra_mex_obj('grad', obj, reshape(x, dimensions)), dimension , 1);


% get starting value gamma
lipschitz_constant= estimate_lipschitz( df,x0 );
gamma0 = (1-beta_safety_value)/lipschitz_constant;
gamma=gamma0;
% gamma = 0.00001;
% %% variables lbgfs
buffer_size=50; % buffer_size

% internal buffers used
alpha=zeros(1,buffer_size);
beta=zeros(1,buffer_size);

s=zeros(dimension,buffer_size); % x_{k+1} - x_{k}
y=zeros(dimension,buffer_size); % df(x_{k+1}) - df(x_{k})
%% iterate
options = struct('GradObj','on','Display','iter','LargeScale','off','HessUpdate','lbfgs', ...
    'InitialHessType','identity','GoalsExactAchieve',1,'GradConstr',false, 'MaxIter', fminlbfgs_iterations);
tic

[exitflag, grad, data, optim] = fminlbfgs_data_init(LM_obj_reg,R_lm0,options);


tau = 1;
bcktrckingSteps = zeros(1,          number_of_steps);


x=x0;% set the starting point, and iterate
for interation_index=1:number_of_steps
    [x_prox_grad_descent,new_gamma] = prox_grad_descent_step( x,gamma,beta_safety_value,proxg,f,df );
    gamma=new_gamma;
    step_prox = (x-x_prox_grad_descent);

    if ~exist('x_lbfgs'), x_prev = x0(:); else, x_prev = data.xInitial; end
    [x_lbfgs,fval_lbfgs,exitflag,output, grad, data] = fminlbfgs_iteration_split1(LM_obj_reg, optim, data, exitflag, 0);

%     R= @(x) (1/gamma)*(x - proxg(x-df(x)*gamma));
%     [s_new,y_new,x_lbfgs] = lbfgs(interation_index,buffer_size,x,R,s,y);
%     s=s_new;y=y_new;
    
    % find the right convex combination trough backtracking
    tau=1;
    max_number_of_steps_backtracking=4;
    for i=1:max_number_of_steps_backtracking
        bcktrckingSteps(interation_index) = bcktrckingSteps(interation_index) + 1 ; 
        d=x_lbfgs(:)-x;
%         step_parameter = 0.01*ceil(interation_index/4);
%         potential_x=x + step_parameter/(norm(step_prox)+norm(d)) * (-(1-tau)*step_prox + tau*d);
        potential_x=x + 1 * (-(1-tau)*step_prox + tau*d);
        sigma=beta_safety_value/(4*gamma);
%         if(FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg  )<= FBE( x,gamma,beta_safety_value,f,df,g,proxg  )-sigma*norm(step_prox/gamma,2)^2)
        if(FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg  ) < 1.2 * FBE( x,gamma,beta_safety_value,f,df,g,proxg))
%             tau=min(1, tau*2)
% %             while(FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg  )<= FBE( x,gamma,beta_safety_value,f,df,g,proxg  )-sigma*norm(step_prox/gamma,2)^2)
%             while(FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg  )<= FBE( x,gamma,beta_safety_value,f,df,g,proxg  ))
%                 tau=min(1, tau*2)
% 
%                 tau=tau*2
%                 potential_x=x + 0.01/(norm(step_prox)+norm(d)) * (-(1-tau)*step_prox + tau*d);
% 
%             end
%             tau=tau/2;
            for get_conditions = 1
            condition1           = FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg ) ;
            condition2_first     = FBE( x,gamma,beta_safety_value,f,df,g,proxg);
            condition2_second    = max(-10, -sigma*norm(step_prox/gamma,2)^2);
            
            condition2_all       = condition2_first + condition2_second;
            condition_bool_first = condition1 > condition2_first;
%             condition_bool_all   = condition1 > condition2_all;
            condition_bool_all   = toc;
            conditionz = [condition1; ...
                         condition2_first; condition2_second; condition2_all; ...
                         condition_bool_first; condition_bool_all];
            end
            
            break; % if this is statified stop right away
        else
            tau=tau/2;
        end
        for get_conditions = 1
            condition1           = FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg)  ;
            condition2_first     = FBE( x,gamma,beta_safety_value,f,df,g,proxg);
            condition2_second    = max(-10, -sigma*norm(step_prox/gamma,2)^2);
            condition2_all       = condition2_first + condition2_second;
            condition_bool_first = condition1 > condition2_first;
%             condition_bool_all   = condition1 > condition2_all;
            condition_bool_all   = toc;
            conditionz = [condition1; ...
                         condition2_first; condition2_second; condition2_all; ...
                         condition_bool_first; condition_bool_all];
            end
            
    end
    condition_array(:, interation_index) = conditionz;

    x                               = potential_x;
    x_steps(:,interation_index)     = x;
    f_evals(interation_index)       = f(x);
    f_lbfgs_evals(interation_index) = fval_lbfgs;
%     f_lbfgs_evals(interation_index) = f(x_lbfgs);
    f_prox_evals(interation_index)  = f(x_prox_grad_descent);
    %     x_actuals(:, interation_index)  = x;
    taus(interation_index)          = tau;
    grad_sizes(:, interation_index) = [norm(x); norm(step_prox); norm(d)];

    updateProxFminLBFGS = 1;
    if updateProxFminLBFGS
            data.xInitial = potential_x;
            data.dir = (data.xInitial - x_prev)/data.alpha;
            data.alpha;
            % UPDATE HESSIAN VIA FminLBFGS    
            [x_lbfgs,fval_lbfgs,exitflag,output, grad, data] = fminlbfgs_iteration_split1(LM_obj_reg, optim, data, exitflag, 1);

    end            

    
end
panoc_t = toc

f_all_evals = [f_evals; f_lbfgs_evals; f_prox_evals; taus]
f_opt = f(Rh);
% f_evals
% f_lbfgs_evals
% bcktrckingSteps
%% plot the convergence rate
figure(1);clf;
subplot(2,1,1)
plot(taus);

% figure(2);clf;
subplot(2,1,2)
plot(f_evals)    
%%
figure;
plot(f_evals, 'ko', 'MarkerSize', 2)
hold on
plot(f_lbfgs_evals, '-.')
plot(f_prox_evals, '--')
plot(condition_array(1,:),'b')
plot(condition_array(2,:),'r')
plot(f_opt*ones(size(f_lbfgs_evals)))

legend('f evals', 'f lbfgs evals', 'f prox evals', 'fbe 1', 'fbe 2', 'optimal f (ground truth)')
%%
figure;
plot(f_evals, 'o', 'MarkerSize', 2, 'MarkerFaceColor','b')
hold on
plot(condition_array(1,:))
plot(condition_array(2,:))
legend('f evals', 'fbe 1', 'fbe 2')
% f(Rh)


% figure(2);clf;
% plot(x_steps(1,2:end)./x_steps(1,1:end-1));