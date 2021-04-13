% COMPARE VARIOUS SOLUTIONS TO IDENTIFICATION PROBLEM
% SPECIFICALLY COMPARE BETWEEN R MATRICES

% From INITIAL SYSTEM
R = ss2r(sys0);
% From IDENT with 'm' solver
Rh = ss2r(sysh_ident);
% From KUNG REALIZATION
Rkung = ss2r(sysh_kung);
% From SLRA_MEX_OBJ('Optimize')
Ropt = info.Rh;

[R ; Rh ; Rkung ; Ropt]

%% RESIDUALS (Distance from Constraint : (R*R' - I = 0)

[C(R) ; C(Rh) ; C(Rkung) ; C(Ropt)]


%% F(R) Values
f_R = slra_mex_obj('func', obj, R);
f_Rh = slra_mex_obj('func', obj, Rh);
f_Rkung = slra_mex_obj('func', obj, Rkung);
f_Ropt = slra_mex_obj('func', obj, Ropt);

[f_R ; f_Rh ; f_Rkung ; f_Ropt]

%% SS's

ss_R = r2ss(R, m, ell);
ss_Rh = r2ss(Rh, m, ell);
ss_Rkung = r2ss(Rkung, m, ell);
ss_Ropt = r2ss(Ropt, m, ell);

sys_comparison(u0, y0, ss_R)
sys_comparison(u0, y0, ss_Rh)
sys_comparison(u0, y0, ss_Rkung)
sys_comparison(u0, y0, ss_Ropt)

