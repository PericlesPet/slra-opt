function [logdata, data_opt, minf_log, x] = ...
    GDescALM(f, df, x0, options)

x         = x0;
maxIter     = options.maxIter; 
gamma       = options.gamma
reg         = options.reg; 
opt         = options.opt;
obj         = options.obj;
R           = options.Ropt; 
extCond     = options.extCond;

clear logdata;
Rini = x;
last_x_elements = 500;
maxIter = ceil(maxIter / last_x_elements) * last_x_elements;
f_log = zeros(maxIter, 1);
minf_log = zeros(ceil(maxIter/100),1);
min_f = inf;

% Gradient Update Parameters
constant_thing = 0.1;   % To avoid division by zero
stepsize_factor = 0.01;  % Initial value of stepsize parameter
stepsize_update = 20;  % Reduce stepsize every X steps
% If the optimal value changes less than (tolPcntChange)% 
% across 100 iterations, stop the descent
tolPcntChange = 1;      
%%
i = 1;

% Stopping criteria (if implemented?)%%
% numberz = 2800;
% (f_log(numberz +100)/f_log(numberz )-1)*100 > 0.1

%%
for i = 1:maxIter

    %% Every 100 iterations print
    if mod(i,100) == 0  
        fprintf('ITERATION %d, minF = %4.4f\n', i, min_f);
        minf_log(i/100) = min_f;
        if (i >= 102 && extCond)        
            fprintf('   PcntChange = %f\n', ...
                abs((minf_log(i/100)/minf_log(i/100-1)-1)*100));
            if abs((minf_log(i/100)/minf_log(i/100-1)-1)*100) < tolPcntChange
                break
            end
        end
    end
    %% Calc f and df
    fval = f(x);
    if i == 2
        min_f = fval ; 
    end
    
    if exist('g') prev_dir = sign(grad); else prev_dir = zeros(size(x)); end
    grad = df(x);
    curr_dir = sign(grad);

    %% Data Logging
    residual_R = x*x' ;
    dir_diff    = sign(R-x) - sign(-grad);
    f_log(i) = fval; 
    id = mod(i-1, last_x_elements)+1;

    logdata(id).f           = fval;
    logdata(id).g           = grad;
    logdata(id).residual_R  = residual_R;
    logdata(id).Rin         = x;
    logdata(id).R_opt       = R;
    logdata(id).R_opt_dist  = R - x;
    logdata(id).dir_diff    = dir_diff;
    logdata(id).R_cond      = cond(x);
    if reg
        logdata(id).f_reg       = f_reg;
        logdata(id).g_reg       = g_reg;
        logdata(id).pcntChange  = norm(gamma*g_reg)/norm(x);
        logdata(id).normG_reg   = norm(gamma*g_reg);
    else
        logdata(id).pcntChange  = norm(gamma*grad)/norm(x);
        logdata(id).normG       = norm(gamma*grad);
    end
    logdata(id).normRin     = norm(x);
    logdata(id).gdDir_diff  = prev_dir - curr_dir;

    %% Descent Step
    stepsize_param = stepsize_factor/ceil(i/stepsize_update);
    %     Rin = Rin - factor*norm(Rin)*gamma * g / (norm(gamma*g) + constant_thing)  ;
    x = x - stepsize_param * gamma * grad / (norm(gamma*grad) + constant_thing)  ;

    %% Kinda search, 
    if reg
        if f_reg <= min_f
            data_opt = logdata(id);
            min_f = f_reg;
        end
    else
        if fval<=min_f
            data_opt = logdata(id);
            min_f = fval;
        end
    end
end
