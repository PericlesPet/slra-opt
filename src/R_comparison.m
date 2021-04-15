% COMPARE VARIOUS SOLUTIONS TO IDENTIFICATION PROBLEM
% SPECIFICALLY COMPARE BETWEEN R MATRICES

% R      : from ACTUAL SYSTEM sys0 (Ground Truth)
% Rh     : from ident Solution
% Rkung  : from Kung Realization
% RmyOpt : from some of my optimization routines??? 
% M_S    : Matrix Structure 
function [R_stats] = R_comparison(Rs, SSs, is_R, C, obj, m, ell, u0, y0, th2R, M_S, psi)

if is_R
    R = Rs.R;
    Rh = Rs.Rh;
    Rkung = Rs.Rkung;
    RmyOpt = Rs.RmyOpt;
else
    ss_R = SSs.sys0;
    ss_Rh = SSs.sysh_ident;
    ss_Rkung = SSs.sysh_kung;
    ss_RmyOpt = SSs.sysh_myOpt;    
    % From INITIAL SYSTEM
    R = ss2r(ss_R);
    % From IDENT with 'm' solver
    Rh = ss2r(ss_Rh);
    % From KUNG REALIZATION
    Rkung = ss2r(ss_Rkung);
    % From SLRA_MEX_OBJ('Optimize')
    RmyOpt = ss2r(ss_RmyOpt);        
end
R_values = [R ; Rh ; Rkung ; RmyOpt]

%% RESIDUALS (Distance from Constraint : (R*R' - I = 0)



th_R      = R2th(R,      M_S, psi, 0);
th_Rh     = R2th(Rh,     M_S, psi, 0);
th_Rkung  = R2th(Rkung,  M_S, psi, 0);
th_RmyOpt = R2th(RmyOpt, M_S, psi, 0);

R_constraints = [C(th_R) ; C(th_Rh) ; C(th_Rkung) ; C(th_RmyOpt)]


%% F(R) Values
f_R = slra_mex_obj('func', obj, R);
f_Rh = slra_mex_obj('func', obj, Rh);
f_Rkung = slra_mex_obj('func', obj, Rkung);
f_Ropt = slra_mex_obj('func', obj, RmyOpt);

R_function_values = [f_R ; f_Rh ; f_Rkung ; f_Ropt]

%% SS's

ss_R = r2ss(R, m, ell);
ss_Rh = r2ss(Rh, m, ell);
ss_Rkung = r2ss(Rkung, m, ell);
ss_Ropt = r2ss(RmyOpt, m, ell);

M_R = sys_comparison(u0, y0, ss_R)
M_Rh = sys_comparison(u0, y0, ss_Rh)
M_Rkung = sys_comparison(u0, y0, ss_Rkung)
M_Ropt = sys_comparison(u0, y0, ss_Ropt)

R_stats = [M_R ; M_Rh ; M_Rkung ; M_Ropt];
