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

%%%% IMPORTANT: RUN ALMTEST FIRST (FOR DEPENDENCIES)

f   = almData.f;
df  = almData.df;
DxL = almData.DxL;
L   = almData.L;
DxLx = @(x) almData.DxL(x, lambda, 10000);
Lx   = @(x) almData.L(x, lambda, 10000);
ce  = almData.ce;
dce = almData.dce;
consGradient = true;
x0 = [ph_ini(:) ; Rini(:)];

outerLoops  = 1;
rho         = 10;
rho_Mult    = 5;

maxItersBase = 15;

maxComplexity = 5;

maxIters     = ceil(maxItersBase * ...
    min((statsTable.complexities(cmplx_iter) / 500),maxComplexity) / 5)*5;

innerLoopsBase  = 8;
innerLoops = ceil(innerLoopsBase * ...
    min((statsTable.complexities(cmplx_iter) / 500),maxComplexity) / 4)*4;




if ~exist('selectFmincons'), selectFmincons = [2:3]; end

for selectAlgo  = selectFmincons

clear fminconData

switch selectAlgo
	case 1 	% ALM function, no constraints (maybe R*R' = I ?)
        fprintf('INITIATE FMINCON for ALM function\n');
		fun     = @(x) objfungrad(x, Lx, DxLx);
%         nonlcon = @(x) [];
        nonlcon = @(x) myConfungrad(x, ce, dce);
        outerLoops = 5;
	case 2  % SLRA function, s.t. R*R' = I
        fprintf('INITIATE FMINCON for SLRA function\n');
		fun     = @(x) LM_obj_reg(reshape(x,size(Rini)));
		ce_rri2 = @(x) stiefConstraint(reshape(x,size(Rini)), 'dist');
		dce_rri2 = @(x)FDGradient(ce_rri2,x,gstep);
		nonlcon = @(x) myConfungrad(x, ce_rri2, dce_rri2);
        x0 = [Rini(:)];
        outerLoops = 1;
	case 4  % norm(p - ph)^2 s.t. constraints,  WITHOUT  Constraint Gradients 
		fun     = @(x) objfungrad(x, f, df);
        dce_empty = @(x)[];
        nonlcon = @(x) myConfungrad(x, ce, dce_empty);
        consGradient = false;
	case 3  % norm(p - ph)^2 s.t. constraints,   WITH  	 Constraint Gradients 
        f   = @(X)1/2*norm(alm_weights(~isinf(alm_weights)).*(X(1:np_w)-wtfdata.p(1:np_w)))^2;
        df  = @(X)[alm_weights(~isinf(alm_weights)).*(X(1:np_w)-wtfdata.p(1:np_w));zeros(R_n*R_m,1)];
        ce  = @(X)([ce_rh0(X),ce_rri(X)]);
        dce = @(X)([dce_rh0(X),dce_rri(X)]);
        fprintf('INITIATE FMINCON for |p - x|^2, s.t. constraints\n');
		fun     = @(x) objfungrad(x, f, df);
        nonlcon = @(x) myConfungrad(x, ce, dce);
        x0 = [ph_ini(:) ; Rini(:)];
        outerLoops = 1;
end

% nonlcon = @(x) confungrad(x, ce_rri, dce_rri);
% nonlcon = @(x) confungrad(x, ce, dce);
dispLvl = 'off';
options = optimoptions('fmincon',...
    'SpecifyObjectiveGradient',true, ...
    'SpecifyConstraintGradient',consGradient, ...
    'Display',dispLvl, ...
    'MaxFunctionEvaluations', 100000, ...
    'MaxIterations', maxIters);
    
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
    DxLx    =   @(x) almData.DxL(x, lambda, rho);
    Lx      =   @(x) almData.L(x, lambda, rho);
    fun     =   @(x) objfungrad(x, Lx, DxLx);
    nonlcon =   @(x) myConfungrad(x, ce, dce);
end

fprintf('   OUTER ITER %d\n', j);

for i = 1:innerLoops
%     fprintf('FMINCON LOOP %d\n', j);
%     fprintf('   inner iter %d\n', i);
%     fprintf('        RHO  %3d\n', rho);
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

if selectAlgo == 1
    fminconData_aLM = fminconData;
    fprintf('Complete: Fmincon ALM in t = %f sec.\n', ...
        fminconData.t_stamps(end));
    fprintf('   SLRA Fmin = %f,\n   Mean Acc = %f %% \n', ...        
        fminconData.f_slra_val(end), mean(fminconData.M0(:,end)));
    fprintf('   Max Mean Acc = %f %% \n', ...        
        max(mean(fminconData.M0(:,:))));
    
elseif selectAlgo == 2
    fminconData_slraVSslra = fminconData;
    fprintf('Complete: Fmincon SLRA in t = %f sec.\n', ...
        fminconData.t_stamps(end));
elseif selectAlgo == 3
    fminconData_gdTrue = fminconData;
    fprintf('Complete: Fmincon |p-x|^2 in t = %f sec.\n', ...
        fminconData.t_stamps(end));
end

end
% beep on
% beep 



