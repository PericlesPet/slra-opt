%% Gradient Stepsize Parameter
gdInput.gamma       = 0.01;
gdInput.Rin         = Rini;
gdInput.maxIter     = 4000;
gdInput.reg         = exist('opt') && isfield(opt, 'method') && strcmp(opt.method, 'reg');
gdInput.reg         = 0;     % Use regularization
gdInput.extCond     = 0; % Use Exit condition to stop iterations
gdInput.obj         = obj; 
gdInput.sysAccuracy = sysAccuracy;
gdInput.Ropt        = R_slramex;

gdAlgo  = 3; % 1 -> simple/reg GD

% f = @(R) slra_mex_obj('func', obj, R);


tic
switch gdAlgo
    case 1   % 1 -> simple/reg GD
        [gdRegData, dataOptId]  = GDesc(gdInput, opt, R_slramex);

        R_gd = reshape(gdRegData.Rin(:,dataOptId), size(Rini));
        gdRegData.Ropt          = R_slramex(:);
        [~, gdRegData.Mopt]     = sysAccuracy(R_slramex);
        gdRegData.fslraOpt      = slra_mex_obj('func', obj, R_slramex);
    case 2   % 2 -> Manifold   GD
        [gdManoptData, dataOptId] = GDescManopt(gdInput, opt, R_slramex);

        R_gdProj = reshape(gdManoptData.Rin(:,dataOptId), size(Rini));
        gdManoptData.Ropt       = R_slramex(:);
        [~, gdManoptData.Mopt]  = sysAccuracy(R_slramex);
        gdManoptData.fslraOpt   = slra_mex_obj('func', obj, R_slramex);
    case 3   % 3 -> Projected  GD
        [gdProjData, dataOptId] = GDescProj(gdInput, opt, R_slramex);

        R_gdProj = reshape(gdProjData.Rin(:,dataOptId), size(Rini));
        gdProjData.Ropt = R_slramex(:);
        [~, gdProjData.Mopt] = sysAccuracy(R_slramex);
        gdProjData.fslraOpt = slra_mex_obj('func', obj, R_slramex);

end

toc


%% Compare all Rs?
% So far:
% R          : slra_mex_obj('optimize') 
% RmyOpt     : optimal 
% Rini       : initial approximation
% sys0       : The Actual random system (ground-truth)
% sysh_ident : result of ident function
% sysh_kung  : result of Kung Realization

sample_divisor1 = 1;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor1),:),y0(1:ceil(length(y0)/sample_divisor1),:), r2ss(R_ident, m_in, ell));
sample_divisor2 = 1;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor2),:),y0(1:ceil(length(y0)/sample_divisor2),:), r2ss(R_gd, m_in, ell));
sample_divisor3 = 1;
sys_comparison(u(1:ceil(length(u)/sample_divisor3),:),y(1:ceil(length(y)/sample_divisor3),:), r2ss(Ropt, m_in, ell));
sample_divisor4 = 1;
sys_comparison(u0(1:ceil(length(u0)/sample_divisor4),:),y0(1:ceil(length(y0)/sample_divisor4),:), r2ss(Rini, m_in, ell));

is_R = 1; 
if is_R  %Use SS format
    Rs.R = ss2r(sys0);
    Rs.Rh = Rh;    
    if ~exist('sysh_kung') Rs.Rkung = Rini; else Rs.Rkung = ss2r(sysh_kung); end
    Rs.RmyOpt = R_gd;    

    SSs = [];
else % Use SS format
    SSs.sys0        = sys0;
    if exist('sysh_ident') SSs.sysh_ident  = sysh_ident; else SSs.sysh_ident = r2ss(Rh); end
    if exist('sysh_kung') SSs.sysh_kung = sysh_kung; else SSs.sysh_kung = r2ss(Rini); end
    SSs.sysh_myOpt  = r2ss(R_gd);

    RRs = [];
end
M_S = phi * (s0 + pext(tts + 1));
R_stats = R_comparison(Rs,SSs, is_R, C, obj, m_in, ell, u0, y0, th2R, M_S, psi);


