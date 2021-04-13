function [f, g] = varproFuncAndGrad(obj, R, mu, reg)


f = slra_mex_obj('func', obj, R)
if reg
    f = f + mu * norm(R * R' - eye(size(R*R')),'fro')^2;
end

if nargout > 1
    g = slra_mex_obj('grad', obj, R)
    if reg
        regularizer_grad =  2 *(Rin * Rin' - eye(size(Rin*Rin')))*Rin;
        g = g + mu * regularizer_grad; 
    end
end


