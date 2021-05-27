options = optimoptions('fmincon',...
    'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);

fun     = @(x) objfungrad(x, f, df);
nonlcon = @(x) confungrad(x, ce, dce);
%%
f   = checkdata.f;
df  = checkdata.df;
DxL = checkdata.DxL;
L   = checkdata.L;
ce  = checkdata.ce;
dce = checkdata.dce;

fcol    = @(x)deal(f(x),df(x));
cecol    = @(x)deal([], ce(x),[], dce(x));

%%

x0;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

%%


function [fval,gradf] = objfungrad(x, f, df)
fval = f(x);
if nargout  > 1
    gradf = df(x);
end
end



function [c,ceq,DC,DCeq] = confungrad(x, ce, dce)
c = [];
% No nonlinear equality constraints
ceq=ce(x);
% Gradient of the constraints:
if nargout > 2
    DC= [];
    DCeq = dce(x);
end
end