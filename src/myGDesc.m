function [logdata, data_opt] = myGDesc(Rin, maxIter, gamma, reg, opt, obj, R)
clear logdata;
last_x_elememts = 500;
maxIter = ceil(maxIter / last_x_elememts) * last_x_elememts;

min_f = inf; 
% Regularizer Parameter
if reg
    mu = opt.g;
end

for i = 1:maxIter
%%
f = slra_mex_obj('func', obj, Rin);
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
g = slra_mex_obj('grad', obj, Rin);
if reg
    regularizer_grad =  2 *(Rin * Rin' - eye(size(Rin*Rin')))*Rin;
    g_reg = g + mu * regularizer_grad; 
end

%%
residual_R = Rin*Rin' ;
dir_diff    = sign(R-Rin) - sign(-g);

if reg    
    Rin = Rin - gamma * g_reg ;
else
    Rin = Rin - gamma * g ;
end

%%
id = mod(i-1, last_x_elememts)+1;

logdata(id).f           = f;
logdata(id).g           = g;
if reg
    logdata(id).f_reg       = f_reg;
    logdata(id).g_reg       = g_reg;
end
logdata(id).residual_R  = residual_R;
logdata(id).Rin         = Rin;
logdata(id).R_opt       = R;
logdata(id).R_opt_dist  = R - Rin;
logdata(id).dir_diff    = dir_diff;
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