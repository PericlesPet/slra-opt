% Simple function to return the value of given ce and dce at x
% Suitable for matlab optimization functions that require "nonlcon" lambda functions
% e.g.
% 		ce = @(x) stiefConstraint(reshape(x,size(Rini)), 'dist'); (Equality Constraint)
% 		dce = @(x)FDGradient(ce,x,gstep);               (Equality Constraint Gradient)
% 		nonlcon = @(x) myConfungrad(x, c, dce);
%       [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

function [c,ceq,DC,DCeq] = myConfungrad(x, ce, dce)
c = [];
% No nonlinear equality constraints
ceq=ce(x);
% Gradient of the constraints:
if nargout > 2
    DC= [];
    DCeq = dce(x);
end
end