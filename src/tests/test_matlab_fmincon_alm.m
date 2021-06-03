%%
f   = checkdata.f;
df  = checkdata.df;
DxL = checkdata.DxL;
L   = checkdata.L;
DxLx = @(x) checkdata.DxL(x, lambda, 10000);
Lx   = @(x) checkdata.L(x, lambda, 10000);
ce  = checkdata.ce;
dce = checkdata.dce;

rho = 10;
consGradient = true;
outerLoops = 1;
x0 = [ph_ini(:) ; Rini(:)];

selectAlgo = 4;
switch selectAlgo
	case 1 	% ALM function, no constraints (maybe R*R' = I ?)
		fun     = @(x) objfungrad(x, Lx, DxLx);
%         nonlcon = @(x) [];
        nonlcon = @(x) myConfungrad(x, ce, dce);        
        outerLoops = 10;
	case 2  % SLRA function, s.t. R*R' = I
		fun     = @(x) LM_obj_reg(reshape(x,size(Rini)));
		ce_rri2 = @(x) stiefConstraint(reshape(x,size(Rini)), 'dist');
		dce_rri2 = @(x)FDGradient(ce_rri2,x,gstep);
		nonlcon = @(x) myConfungrad(x, ce_rri2, dce_rri2);
        x0 = [Rini(:)];
	case 3  % norm(p - ph)^2 s.t. constraints,  WITHOUT  Constraint Gradients 
		fun     = @(x) objfungrad(x, f, df);
        dce_empty = @(x)[];
        nonlcon = @(x) myConfungrad(x, ce, dce_empty);
        consGradient = false;
	case 4  % norm(p - ph)^2 s.t. constraints,   WITH  	 Constraint Gradients 
		fun     = @(x) objfungrad(x, f, df);
        nonlcon = @(x) myConfungrad(x, ce, dce);        

end


% ce_rri
% dce_rri


% nonlcon = @(x) confungrad(x, ce_rri, dce_rri);
% nonlcon = @(x) confungrad(x, ce, dce);
maxIters = 10;
options = optimoptions('fmincon',...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',consGradient, ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 100000, ...
    'MaxIterations', maxIters);
    

%%
clear fminconData;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x0_temp = x0;

index = 1;
rho = 10;

for j = 1:outerLoops
if selectAlgo == 1
    DxLx    =   @(x) checkdata.DxL(x, lambda, rho);
    Lx      =   @(x) checkdata.L(x, lambda, rho);
    fun     =   @(x) objfungrad(x, Lx, DxLx);
    nonlcon =   @(x) myConfungrad(x, ce, dce);
end

for i = 1:500
%     fprintf('FMINCON LOOP %d\n', j);
    fprintf('        ITER %d\n', i);
    fprintf('        RHO  %3d\n', rho);
    tic
%     if i==1, options.MaxIterations=1; else options.MaxIterations = maxIters;end
    
    [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

    if exist('fminconData'), t_stamp = fminconData.t_stamps(end) + toc; ...
    else, t_stamp = toc; end


    if length(x) == numberofelements(slradata.Rini)
        R_fmincon = reshape(x, size(slradata.Rini));       
    else
        R_fmincon = reshape(x(slradata.npExt+1:end), size(slradata.Rini));    
    end
    
    [~, M0] = compare(iddata(slradata.y0, slradata.u0), ...
        idss(r2ss(R_fmincon, slradata.m_in, slradata.ell))); 

    fminconData.Xs(:, index)        = x;
    fminconData.fvals(index)        = fval; 
    fminconData.M0(index)           = mean(M0);
    fminconData.f_slra_val(index)   = slra_mex_obj('func', obj, R_fmincon);
%     fminconData.CE(:,index)         = ce(x);    
%     fminconData.CE(:,index)         = ce(x0);
    [~, fminconData.CE(:,index), ~, ~] = nonlcon(x);
    fminconData.t_stamps(index)     = t_stamp;
    x0 = x;
    index = index + 1;    
end

rho = rho * 5;
end

x0 = x0_temp;

plot(fminconData.f_slra_val)
min(fminconData.f_slra_val)
max(fminconData.M0)
beep on
beep 
%%
figure
subplot(2,2,1)
plot(fminconData.t_stamps, fminconData.f_slra_val)
title('SLRA M(R)')

subplot(2,2,2)
semilogy(fminconData.t_stamps, fminconData.CE(1,:))
% hold on 
% semilogy(fminconData.CE(2,:))
title('Constraints')

subplot(2,2,3)
plot(fminconData.t_stamps, fminconData.fvals)
title('Fvals')

subplot(2,2,4)
plot(fminconData.t_stamps, fminconData.M0)
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


function [fval,gradf] = objfungrad(x, f, df)
fval = f(x);
if nargout  > 1
    gradf = df(x);
end
end



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