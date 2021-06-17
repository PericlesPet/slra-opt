%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This requires having already ran ALM (almTest, almSearch ...)
% Check slra M(R), L(x,R,Y), f(x) = norm(p-x)^2 THRU matlab's fmincon (nonlinear optimization)
% THIS SCRIPT PERFORMS OPTIMIZATION USING MATLAB'S FMINCON FUNCTION 
% FOR VARIOUS SETTINGS:
% 1: Optimize the ALM function, given the stiefel constraint R*R' = I 
%       L(x) = @(x) L(x,lambda, c) = ...
% 
% 2: Optimize the SLRA Function, given the stiefel constraint R*R' = I
%       M(R) = min( norm(p-phat)^2 ) over phat,   s.t. R*R' = I
% 
% 3: Optimize the initial function before VARPRO
%       min( norm(p-phat)^2 ) over phat, R,   s.t. R*Hank(phat) = 0
%                                                  R*R'         = I
% 
% 4: Same as 3 but provide the Constraint Gradients 
%    and don't let Matlab calculate them
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
selectAlgo  = 2;


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

    % GATHER DATA    
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
% ----> dataVisualize 2ND PART



%% CHECK SLRA ACCURACY
for i = 1:size(info.RhK, 3)
    [~, M_testslra(:,i)] = sysAccuracy(info.RhK(:,:,i));    
end
subplot(2,1,1)
plot(mean(M_testslra))
subplot(2,1,2)
plot(info.iterinfo(1,:), info.iterinfo(2,:), 'r--')

max(mean(M_testslra))






