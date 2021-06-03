%% Check slra M(R), L(x,R,Y), f(x) = norm(p-x)^2 THRU matlab's fmincon (nonlinear optimization)
% This requires having already ran 
f   = checkdata.f;
df  = checkdata.df;
DxL = checkdata.DxL;
L   = checkdata.L;
DxLx = @(x) checkdata.DxL(x, lambda, 10000);
Lx   = @(x) checkdata.L(x, lambda, 10000);
ce  = checkdata.ce;
dce = checkdata.dce;
consGradient = true;
x0 = [ph_ini(:) ; Rini(:)];

outerLoops  = 1;
rho         = 10;
rho_Mult    = 5;

maxIters    = 8;
innerLoops  = 10;
selectAlgo  = 3;


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
options = optimoptions('fmincon',...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',consGradient, ...
    'Display','iter', ...
    'MaxFunctionEvaluations', 100000, ...
    'MaxIterations', maxIters);
    

%
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

for i = 1:innerLoops
%     fprintf('FMINCON LOOP %d\n', j);
    fprintf('        ITER %d\n', i);
    fprintf('        RHO  %3d\n', rho);
    tic
    if i==1, options.MaxIterations=1; else options.MaxIterations = maxIters;end
    
    [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

    if exist('fminconData'), t_stamp = fminconData.t_stamps(end) + toc; ...
    else, t_stamp = toc; end


    if length(x) == numel(slradata.Rini)
        R_fmincon = reshape(x, size(slradata.Rini));       
    else
        R_fmincon = reshape(x(slradata.npExt+1:end), size(slradata.Rini));    
    end
    
    [~, M0] = compare(iddata(slradata.y0, slradata.u0), ...
        idss(r2ss(R_fmincon, slradata.m_in, slradata.ell))); 

    fminconData.Xs(:, index)            = x;
    fminconData.fvals(index)            = fval; 
%     fminconData.M0(index)               = mean(M0);
    fminconData.M0(:,index)               = M0;
    fminconData.f_slra_val(index)       = slra_mex_obj('func', obj, R_fmincon);
    %     fminconData.CE(:,index)           = ce(x);    
    %     fminconData.CE(:,index)           = ce(x0);
    [~, fminconData.CE(:,index), ~, ~]  = nonlcon(x);
    fminconData.t_stamps(index)         = t_stamp;
    x0 = x;
    index = index + 1;    
end

rho = rho * rho_Mult;
end

x0 = x0_temp;

plot(fminconData.f_slra_val)
min(fminconData.f_slra_val)
max(mean((fminconData.M0)))
beep on
beep 
%% VISUALIZE FVALS, SLRA, CONSTRAINTS, AND ACCURACY. ALSO COMPARE TO OPTIMAL SLRA_MEX VALUES

% fminconData = fminconData_aLM;
% fminconData = fminconData_gdTrue;
% fminconData = fminconData_gdfalse;
fminconData = fminconData_slraVSslra
use_iter = 0;
if ~isfield(fminconData, 't_stamps'), ...
        fminconData.t_stamps = 1:length(fminconData.f_slra_val); use_iter = 1;end


set(0,'DefaultLegendAutoUpdate','off')
figure
subplot(2,2,1)
plot(fminconData.t_stamps, fminconData.f_slra_val)
title('SLRA M(R)')
hold on
if ~use_iter, plot(info.iterinfo(1,:), info.iterinfo(2,:), 'r--'), else ...
    plot(1:length(info.iterinfo(1,:)), info.iterinfo(2,:), 'r--'), end
line([0 fminconData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')
legend('fmincon', 'slra mex', 'optimum point')

subplot(2,2,2)
for ii = 1:size(fminconData.CE, 1)    
    semilogy(fminconData.t_stamps, fminconData.CE(ii,:))
    hold on
end
% hold on
% semilogy(fminconData.CE(2,:))
title('Constraints')
if size(fminconData.CE, 1) == 1
    legend('RR^T - I_N = 0 Const.')
else
    legend('R*Hank(p_{hat}) = 0 Const.', 'RR^T - I_N = 0 Const.')
end

subplot(2,2,3)
plot(fminconData.t_stamps, fminconData.fvals)
title('Fvals')

subplot(2,2,4)
plot(fminconData.t_stamps, max(fminconData.M0, -100), 'b')
% ylim([-200 100])
title('Accuracy (%)')
if size(fminconData.M0, 1) == 2
    line([0 fminconData.t_stamps(end)], ...
        [M_slra(1) M_slra(1)], 'Color', 'r')
    line([0 fminconData.t_stamps(end)], ...
        [M_slra(2) M_slra(2)], 'Color', 'r')
    legend('accuracy for y1', 'accuracy for y2', 'slra best accuracy for y1','slra best accuracy for y2')
else
    line([0 fminconData.t_stamps(end)], ...
        [mean(M_slra) mean(M_slra)], 'Color', 'r')
    legend('accuracy for fmincon', 'slra best accuracy')
end

% subplot(2,2,3)
% plot(fminconData.f_slra_val)

%% VISUALIZE DIFFERENCE BETWEEN FVALS, SLRA, AND CONSTRAINTS (PROBABLY USELESS)

figure
%%%%% SLRA %%%%%
subplot(2,2,1)
plot(fminconData.t_stamps, fminconData.f_slra_val)
title('SLRA M(R)')
hold on
if ~use_iter, plot(info.iterinfo(1,:), info.iterinfo(2,:), 'r--'), else ...
    plot(1:length(info.iterinfo(1,:)), info.iterinfo(2,:), 'r--'), end
line([0 fminconData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')
legend('fmincon', 'slra mex', 'optimum point')

%%%%% FVALS %%%%%
subplot(2,2,3)
plot(fminconData.t_stamps, 2*fminconData.fvals)
title('Fvals')

%%%%% SLRA - FVALS %%%%%
subplot(2,2,2)
plot(fminconData.t_stamps, fminconData.f_slra_val-2*fminconData.fvals)
title('SLRA - Fvals')

%%%%% CONSTRAINTS%%%%%
subplot(2,2,4)
for ii = 1:size(fminconData.CE, 1)    
    semilogy(fminconData.t_stamps, fminconData.CE(ii,:))
    hold on
end
% hold on
% semilogy(fminconData.CE(2,:))
title('Constraints')
if size(fminconData.CE, 1) == 1
    legend('RR^T - I_N = 0 Const.')
else
    legend('R*Hank(p_{hat}) = 0 Const.', 'RR^T - I_N = 0 Const.')
end




%% CHECK SLRA ACCURACY
for i = 1:size(info.RhK, 3)
    [~, M_testslra(:,i)] = sysAccuracy(info.RhK(:,:,i));    
end
subplot(2,1,1)
plot(mean(M_testslra))
subplot(2,1,2)
plot(info.iterinfo(1,:), info.iterinfo(2,:), 'r--')

max(mean(M_testslra))



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