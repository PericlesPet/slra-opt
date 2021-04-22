
g = slra_mex_obj('grad', obj, Rini)
g_proj = stiefel_proj(Rini, g)
g * g'
g_proj * g_proj'
% stiefel_proj(Rini, g) * stiefel_proj(Rini, g)'
%%
g_proj2 = stiefel_proj(Rini, g_proj2)
g_proj2 * g_proj2'
%%
g_proj_test = projection(Rini, g)
g_proj_test * g_proj_test'

%%
stepsize_param = 0.01 ; 
%%
Rnew  = Rini - stepsize_param/10000000 * g_proj 
Rnew2 = Rini - stepsize_param * g_proj_test / (norm(g_proj_test))  
Rnew3 = Rini - stepsize_param * g / (norm(g))  
%%
Rnew * Rnew' 
Rnew2 * Rnew2' 
Rnew3 * Rnew3' 
Rini * Rini'


%%
f = slra_mex_obj('func', obj, Rnew)
f = slra_mex_obj('func', obj, Rnew2)
f = slra_mex_obj('func', obj, Rnew3)