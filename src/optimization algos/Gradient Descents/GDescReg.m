function [logdata, data_opt, f_log] = GDescReg(Rin, maxIter, gamma, reg, opt, obj, R)
%%
clear logdata;

Rini = Rin;
last_x_elements = 500;
maxIter = ceil(maxIter / last_x_elements) * last_x_elements;
f_log = zeros(maxIter, 1);
min_f = inf; 
% Regularizer Parameter
if reg
    mu = opt.g;
end

% Gradient Update Parameters
constant_thing = 0.1;   % To avoid division by zero
stepsize_factor = 0.1;  % Initial value of stepsize parameter
stepsize_update = 100;  % Reduce stepsize every X steps

%%
i = 1;

%%
for i = 1:maxIter
%%
if mod(i,100) == 0  
    i
    min_f
end
% i = i + 1;
f = slra_mex_obj('func', obj, Rin);
% Rin;
if reg
    regularizer = norm(Rin * Rin' - eye(size(Rin*Rin')),'fro')^2;
    f_reg = f + mu * regularizer;
end

if i == 2
    if reg
        min_f = f_reg;
    else
        min_f = f ; 
    end
end


%%
if exist('g') prev_dir = sign(g); else prev_dir = zeros(size(Rin)); end
g = slra_mex_obj('grad', obj, Rin);
curr_dir = sign(g);
if reg
    regularizer_grad =  2 *(Rin * Rin' - eye(size(Rin*Rin')))*Rin;
    g_reg = g + mu * regularizer_grad ;
end

%%
residual_R = Rin*Rin' ;
dir_diff    = sign(R-Rin) - sign(-g);

%%
f_log(i) = f; 
id = mod(i-1, last_x_elements)+1;

logdata(id).f           = f;
logdata(id).g           = g;
logdata(id).residual_R  = residual_R;
logdata(id).Rin         = Rin;
logdata(id).R_opt       = R;
logdata(id).R_opt_dist  = R - Rin;
logdata(id).dir_diff    = dir_diff;
logdata(id).R_cond      = cond(Rin);
if reg
    logdata(id).f_reg       = f_reg;
    logdata(id).g_reg       = g_reg;
    logdata(id).pcntChange  = norm(gamma*g_reg)/norm(Rin);
    logdata(id).normG_reg   = norm(gamma*g_reg);
else
    logdata(id).pcntChange  = norm(gamma*g)/norm(Rin);
    logdata(id).normG       = norm(gamma*g);
end
logdata(id).normRin     = norm(Rin);
logdata(id).gdDir_diff  = prev_dir - curr_dir;
%%
%%
stepsize_param = stepsize_factor/ceil(i/stepsize_update);

if reg    
%     Rin = Rin - factor*(norm(Rin)+constant_thing)*(gamma * g_reg) / (norm(gamma*g_reg)+constant_thing) ;
    Rin = Rin - stepsize_param * (gamma * g_reg) / (norm(gamma*g_reg)+constant_thing) ;
else
%     Rin = Rin - factor*norm(Rin)*gamma * g / (norm(gamma*g) + constant_thing)  ;
    Rin = Rin - stepsize_param * gamma * g / (norm(gamma*g) + constant_thing)  ;
end

%%
if reg
    if f_reg <= min_f
        data_opt = logdata(id);
        min_f = f_reg;
    end
else
    if f<=min_f
        data_opt = logdata(id);
        min_f = f;
    end
end
end