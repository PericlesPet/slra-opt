%% Gradient Stepsize Parameter
gdInput.gamma       = 0.01;
gdInput.Rin         = Rini;
gdInput.maxIter     = ceil(4000 * min((statsTable.complexities(cmplx_iter) / 500)^1.5, maxComplexity) / 1000)*1000;
gdInput.reg         = exist('opt') && isfield(opt, 'method') && strcmp(opt.method, 'reg');
gdInput.reg         = 0;     % Use regularization
gdInput.extCond     = 0; % Use Exit condition to stop iterations
gdInput.obj         = obj; 
gdInput.sysAccuracy = sysAccuracy;
gdInput.Ropt        = R_slramex;
if ~exist('dspLvlgd'),  dspLvlgd  = 'iter';end
if ~exist('dspFreqgd'), dspFreqgd = 100;end
gdInput.dspLvl = dspLvlgd;
gdInput.dspFreq = dspFreqgd;

if ~exist('selectGDs'), selectGDs = 1:4; end

for gdAlgo  = selectGDs % 1 -> simple/reg GD

tic
switch gdAlgo
    case 1   % 1 -> simple
        gdInput.reg         = 0; 
        [gdData, dataOptId]  = GDesc(gdInput, opt, R_slramex);
        t_gd = toc;
        R_gd = reshape(gdData.Rin(:,dataOptId), size(Rini));
        gdData.Ropt          = R_slramex(:);
        [~, gdData.Mopt]     = sysAccuracy(R_slramex);
        gdData.fslraOpt      = slra_mex_obj('func', obj, R_slramex);
        fprintf('Complete: Simple GD, t = %f\n', t_gd);
        fprintf('   Optimal F : %f \n', gdData.fslraOpt);

    case 2   % 2 -> Manifold   GD
        [gdManoptData, dataOptId] = GDescManopt(gdInput, opt, R_slramex);
        t_gd = toc;
        R_gdManopt = reshape(gdManoptData.Rin(:,dataOptId), size(Rini));
        gdManoptData.Ropt       = R_slramex(:);
        [~, gdManoptData.Mopt]  = sysAccuracy(R_slramex);
        gdManoptData.fslraOpt   = slra_mex_obj('func', obj, R_slramex);
        fprintf('Complete: Manopt-Like GD t = %f\n', t_gd);
        fprintf('   Optimal F : %f \n', gdManoptData.fslraOpt);
        
    case 3   % 3 -> Projected  GD
        [gdProjData, dataOptId] = GDescProj(gdInput, opt, R_slramex);       
        t_gd = toc;
        R_gdProj = reshape(gdProjData.Rin(:,dataOptId), size(Rini));
        gdProjData.Ropt = R_slramex(:);
        [~, gdProjData.Mopt] = sysAccuracy(R_slramex);
        gdProjData.fslraOpt = slra_mex_obj('func', obj, R_slramex);
        fprintf('Complete: Projected GD t = %f\n', t_gd);
        fprintf('   Optimal F : %f \n', gdProjData.fslraOpt);
        
    case 4   % 4 -> reg GD
        gdInput.reg         = 1; 
        [gdRegData, dataOptId]  = GDesc(gdInput, opt, R_slramex);
        t_gd = toc;
        R_gdReg = reshape(gdRegData.Rin(:,dataOptId), size(Rini));
        gdRegData.Ropt          = R_slramex(:);
        [~, gdRegData.Mopt]     = sysAccuracy(R_slramex);
        gdRegData.fslraOpt      = slra_mex_obj('func', obj, R_slramex);   
        fprintf('Complete: Regularized GD t = %f\n', t_gd);
        fprintf('   Optimal F : %f \n', gdRegData.fslraOpt);
end
end