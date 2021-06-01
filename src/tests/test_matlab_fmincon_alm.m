%%
f   = checkdata.f;
df  = checkdata.df;
DxL = checkdata.DxL;
L   = checkdata.L;
DxLx = @(x) checkdata.DxL(x, lambda, 10);
Lx   = @(x) checkdata.L(x, lambda, 10);
ce  = checkdata.ce;
dce = checkdata.dce;

rho = 10;
lambda;

selectAlgo = 1;

switch selectAlgo
	case 1 	% ALM function, no constraints (maybe R*R' = I ?)
		fun     = @(x) objfungrad(x, Lx, DxLx);

	case 2  % SLRA function, s.t. R*R' = I
		fun     = @(x) LM_obj_reg(reshape(x,size(Rini)));
		
		ce_rri2 = @(x) stiefConstraint(reshape(x,size(Rini)), 'dist');
		dce_rri2 = @(x)FDGradient(ce_rri2,x,gstep);
		nonlcon = @(x) confungrad(x, ce_rri2, dce_rri2);

	case 3  % norm(p - ph)^2 s.t. constraints,  WITHOUT  Constraint Gradients 
		fun     = @(x) objfungrad(x, f, df);

	case 4  % norm(p - ph)^2 s.t. constraints,   WITH  	 Constraint Gradients 
		fun     = @(x) objfungrad(x, f, df);

end
% ce_rri
% dce_rri


% nonlcon = @(x) confungrad(x, ce_rri, dce_rri);
% nonlcon = @(x) confungrad(x, ce, dce);

options = optimoptions('fmincon',...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',true, ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 100000, ...
    'MaxIterations', 1);
    

%%
clear fminconData;
% x0 = [ph_ini(:) ; Rini(:)];
x0 = [Rini(:)];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x0_temp = x0;

index = 1;
% rho = 10;
% for j = 1:10
% 
% DxLx = @(x) checkdata.DxL(x, lambda, rho);
% Lx   = @(x) checkdata.L(x, lambda, rho);
% fun     = @(x) objfungrad(x, Lx, DxLx);
% nonlcon = @(x) confungrad(x, ce, dce);

for i = 1:1000
%     fprintf('FMINCON LOOP %d\n', j);
    fprintf('        ITER %d\n', i);
    fprintf('        RHO  %3d\n', rho);
    [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    
%     R_fmincon = reshape(x(slradata.np+1:end), size(slradata.Rini));
    R_fmincon = reshape(x, size(slradata.Rini));
    [~, M0] = compare(iddata(slradata.y0, slradata.u0), ...
        idss(r2ss(R_fmincon, slradata.m_in, slradata.ell))); 

    fminconData.Xs(:, index)        = x;
    fminconData.fvals(index)        = fval; 
    fminconData.M0(index)           = mean(M0);
    fminconData.f_slra_val(index)   = slra_mex_obj('func', obj, R_fmincon);
%     fminconData.CE(:,index)         = ce(x);    
    fminconData.CE(:,index)         = ce_rri2(x);    
    x0 = x;
    index = index + 1;    
end

% rho = rho * 5;

% end
x0 = x0_temp;

plot(fminconData.f_slra_val)
min(fminconData.f_slra_val)
max(fminconData.M0)
beep on
beep 
%%
figure
subplot(2,2,1)
plot(fminconData.f_slra_val)
title('SLRA M(R)')

subplot(2,2,2)
semilogy(fminconData.CE(1,:))
% hold on 
% semilogy(fminconData.CE(2,:))
title('Constraints')

subplot(2,2,3)
plot(fminconData.fvals)
title('Fvals')

subplot(2,2,4)
plot(fminconData.M0)
% ylim([-200 100])
title('Accuracy (%)')
% subplot(2,2,3)
% plot(fminconData.f_slra_val)
% 


%%
R_fmincon = reshape(x(slradata.np+1:end), size(slradata.Rini));
[~, M0] = compare(iddata(slradata.y0, slradata.u0), ...
    idss(r2ss(R_fmincon, slradata.m_in, slradata.ell))) 
fminconData.M0(i) = mean(M0);
fminconData.f_slra_val(i) = slra_mex_obj('func', obj, R_fmincon);
fminconData.CE(:,i) = ce(x);
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