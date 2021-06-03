% Simple function to return the value of given f and df at x
% Suitable for matlab optimization functions that require "fun" lambda functions
% e.g.
%        DxLx = @(x) checkdata.DxL(x, lambda, 10000);
%        Lx   = @(x) checkdata.L(x, lambda, 10000);
% -----> fun     =   @(x) objfungrad(x, Lx, DxLx);  <-----------------
%        [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
function [fval,gradf] = objfungrad(x, f, df)
    fval = f(x);
    if nargout  > 1
        gradf = df(x);
    end
end
