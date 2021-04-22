function [f, g] = varproFuncAndGrad(obj, R, mu, reg, th_format)


f = slra_mex_obj('func', obj, R);
if reg
    f = f + mu * norm(R * R' - eye(size(R*R')),'fro')^2;
end

if nargout > 1
    g = slra_mex_obj('grad', obj, R);
    if reg
        regularizer_grad =  2 *(R * R' - eye(size(R*R')))*R;
        g = g + mu * regularizer_grad; 
    end
    
    if nargin == 5        
        if th_format.active
            PhiS_mat = th_format.PhiS_mat;
            psi = th_format.psi;        
            th_g = R2th(g, PhiS_mat, psi, 0);
            g = th_g;
        end
    end    
end


