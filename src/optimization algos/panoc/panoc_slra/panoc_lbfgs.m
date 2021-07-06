% PANOC algorithm using the lbgfs and proximal gradient methods
% include the Matlab folder with all its subfolders to your path
if isCloseAll == 1
    close all
    clc
end
%maxComplexity       = 15;
number_of_stepsBase = 600;
number_of_steps     = ceil(number_of_stepsBase * ...
    min((statsTable.complexities(cmplx_iter) / 500)^1.5,maxComplexity) / 100)*100;

printFreq       = 100;
dimensions      = size(Rini);
dimension       = dimensions(1) * dimensions(2);
x_steps         = zeros(dimension,  number_of_steps);
grad_sizes      = zeros(3,          number_of_steps);
f_evals         = zeros(1,          number_of_steps);
f_evals_lbfgs   = zeros(1,          number_of_steps);
f_evals_prox    = zeros(1,          number_of_steps);
bcktrckingSteps = zeros(1,          number_of_steps);
taus            = zeros(1,number_of_steps);
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
alpha   =zeros(1,buffer_size);
beta    =zeros(1,buffer_size);
s_lbfgs = zeros(dimension,buffer_size); % x_{k+1} - x_{k}
y_lbfgs = zeros(dimension,buffer_size); % df(x_{k+1}) - df(x_{k})

% iterate
bcktrckingSteps = zeros(1, number_of_steps);
tau             = 1;

index           = 1;
clear panocLbfgsData
nonlcon = @(x) myConfungrad(x, constraints.ce_rri2, constraints.dce_rri2);

x=x0;% set the starting point, and iterate
for interation_index=1:number_of_steps
    tic
    
    [x_prox_grad_descent,new_gamma] = prox_grad_descent_step( x,gamma,beta_safety_value,proxg,f,df );
    gamma=new_gamma;
    step_prox = (x-x_prox_grad_descent);
    
    R= @(x) (1/gamma)*(x - proxg(x-df(x)*gamma));
    [s_new,y_new,x_lbfgs] = lbfgs(interation_index,buffer_size,x,R,s_lbfgs,y_lbfgs);
    s_lbfgs=s_new;y_lbfgs=y_new;
    
    % find the right convex combination trough backtracking
    tau=0;
    max_number_of_steps_backtracking=1;
    for i=1:max_number_of_steps_backtracking
        bcktrckingSteps(interation_index) = bcktrckingSteps(interation_index) + 1 ; 
        d=x_lbfgs-x;
        step_parameter = 0.01*ceil(interation_index/4);
%         potential_x=x + step_parameter/(norm(step_prox)+norm(d)) * (-(1-tau)*step_prox + tau*d);
        potential_x=x + 1 * (-(1-tau)*step_prox + tau*d);
        sigma=beta_safety_value/(4*gamma);
        if(FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg  )<= FBE( x,gamma,beta_safety_value,f,df,g,proxg  )-sigma*norm(step_prox/gamma,2)^2)
%         if(FBE( potential_x,gamma,beta_safety_value,f,df,g,proxg  )<= FBE( x,gamma,beta_safety_value,f,df,g,proxg))
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
            break; % if this is statified stop right away
        else
            tau=tau/2;
        end
    end
    
    
    x                               = potential_x;
%     x_steps(:,interation_index)     = x;
%     f_evals(interation_index)       = f(x);
%     f_evals_lbfgs(interation_index) = f(x_lbfgs);
%     f_evals_prox(interation_index)  = f(x_prox_grad_descent);
    %     x_actuals(:, interation_index)  = x;
    taus(interation_index)          = tau;
    grad_sizes(:, interation_index) = [norm(x); norm(step_prox); norm(d)];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LOG DATA
    if exist('panocLbfgsData'), t_stamp = panocLbfgsData.t_stamps(end) + toc; ...
    else, t_stamp = toc; end
    R_panoclbfgs = reshape(x, size(slradata.Rini));     
    [~, M0] = compare(iddata(slradata.y0, slradata.u0), ...
        idss(r2ss(R_panoclbfgs, slradata.m_in, slradata.ell))); 

    panocLbfgsData.Xs(:, index)               = potential_x(:);
    panocLbfgsData.Xs_lbfgs(:, index)         = x_lbfgs(:);
    panocLbfgsData.Xs_proxgd(:, index)        = x_prox_grad_descent(:);

    panocLbfgsData.fvals(index)               = f(x); 
    panocLbfgsData.fvals_lbfgs(index)         = f(x_lbfgs); 
    panocLbfgsData.fvals_prox(index)          = f(x_prox_grad_descent); 
    panocLbfgsData.f_slra_val(index)          = slra_mex_obj('func', obj, R_panoclbfgs);

    [~, panocLbfgsData.CE(:,index), ~, ~]     = nonlcon(x(:));
    panocLbfgsData.M0(:,index)                = M0;
    panocLbfgsData.t_stamps(index)            = t_stamp;
    %     fminconData.M0(index)               = mean(M0);
    %         fminconData.CE(:,index)           = ce(x);    
    %         fminconData.CE(:,index)           = ce(x0); 

    index = index + 1;         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if mod(interation_index, printFreq) == 0
        fprintf('f(%d) = %f  -  t = %f  -  Accuracy = %f\n', ...
            interation_index, panocLbfgsData.fvals(interation_index),t_stamp, mean(M0)); end


end

if plotPANOC 
    figure
    subplot(2,1,1)
    plot(panocLbfgsData.t_stamps,panocLbfgsData.fvals, 'k:', ...
        'Marker','*', 'MarkerSize', 3, 'MarkerIndices', 1:10:length(panocLbfgsData.t_stamps))
    hold on
    plot(panocLbfgsData.t_stamps,panocLbfgsData.fvals_lbfgs, 'r--')
    plot(panocLbfgsData.t_stamps,panocLbfgsData.fvals_prox, 'b--')
    legend('fvals', 'fvals_{lbfgs}', 'fvals_{prox}')
    ylim([0.9*f(R_true) max(panocLbfgsData.fvals)])
    title('F Evaluations')


    subplot(2,1,2)
    if ~isAccSemilog
        plot(panocLbfgsData.t_stamps,max(mean(panocLbfgsData.M0),0))
    else
        semilogy(panocLbfgsData.t_stamps, ...
        max(mean(panocLbfgsData.M0(:,:)), ...
            1./abs(mean(panocLbfgsData.M0(:,:)))))

        ylim([0 100])
    %     semilogy(panocLbfgsData.t_stamps, mean(panocLbfgsData.M0))
    end
    title('Mean Accuracy')
    suptitle('PANOC_{LBFGS} Figures')
end

fprintf("Time Elapsed on PANOC_LBFGS : %.3f\n", t_stamp);
x_panoclbfgs_final = x_steps(:,end);
R_panoclbfgs       = reshape(x_panoclbfgs_final, size(Rini));