% CALCULATE GRADIENT OF FUNCTION WITH FORWARD DIFFERENCES
function grad = FDGradient(func, x, gstep)

grad=zeros(length(x),1);
fval = func(x);
if nargin < 3 
    gstep=1/1e6;
end

% for i = 1:length(fval)
%     if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
%     if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end

for j=1:length(x)    
    X_temp=x; X_temp(j)=X_temp(j)+gstep;
    [fval_g]=func(X_temp);
    grad(j)=(fval_g-fval)/gstep;
end

end
% end
