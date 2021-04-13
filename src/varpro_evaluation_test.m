% Perform up to "Data Generation" of slra_sw_example

[w, s, r, opt, q, N, T] = ident_preprocessing(w, m_in, ell, opt_oe);
%%
p = w2p(w);
R = ss2r(sys0); 

%%
% tts  - structure specification S(p) = phi * (s0 + p(tts)) 
[tts, p, r, w, Rini, phi ,psi, opt, th2R, C, s0] = slra_to_slraext(p, s, r, opt);

%%
% pext = [0; p];
% prob.x0 = R2th(Rini, phi * (s0 + pext(tts + 1)), psi, opt.R0); 
 
  
%% SLRA HANDLERS

for missing_data = 1
    
    Im = find(isnan(p));
    if ~isempty(Im)
      if ~isfield(s, 'w'), s.w = ones(size(p)); end 
      q = length(s.m); if exist('p'), 
                         np = length(p); 
                       else
                         np = sum(s.m) * length(s.n) + length(s.m) * sum(s.n) ...
                                                     - length(s.m) * length(s.n);
                       end, if ~isfield(s, 'n'), s.n = np - sum(s.m) + 1; end
      N = length(s.n); n  = sum(s.n);
      if length(s.w(:)) == q || all(size(s.w) == [q N])
        % convert q x 1 s.w to q x N
        if isvector(s.w), s.w = s.w(:); s.w = s.w(:, ones(1, N)); end
        % convert q x N s.w to np x 1
        w = [];
        for j = 1:N
          for i = 1:q
            wij = s.w(i, j) * ones(s.m(i) + s.n(j) - 1, 1); w = [w; wij];
          end
        end
        s.w = w;
      end
      p(Im) = 0; s.w(Im) = 0;
    end
end

%%%%% IF SOLVER = C %%%%%%
if isfield(s, 'w'), s.w(find(s.w == 0)) = opt.tol_m; end
opt = rmfield(opt, 'solver');

%% SLRA_MEX_OBJ 

obj = slra_mex_obj('new', p, s, r);
%%
[ph, info] = slra_mex_obj('optimize', obj, opt);

%%
R = info.Rh;
% Rin = ones(size(R));
Rin = Rini;
%%
% Gradient Stepsize Parameter
gamma = 0.0005;

% Regularizer Parameter
opt.g = norm(p(Ig)) ^ 2;
mu = opt.g;

%% 
clear logdata;
for ii = 1:100
%%
f = slra_mex_obj('func', obj, Rin);
% f = slra_mex_obj('func', obj, R);

regulariz0r = norm(Rin * Rin' - eye(size(Rin*Rin')),'fro')^2;
f_reg = f + mu * regulariz0r;

% prob.objective = @(th) slra_mex_obj('func', obj, th2R(th)) ...
%                          + mu * norm(C(th), 'fro') ^ 2;
                     
                     
%%
g = slra_mex_obj('grad', obj, Rin);
regulariz0r_grad =  2 *(Rin * Rin' - eye(size(Rin*Rin')))*Rin;
g_reg = g + mu * regulariz0r_grad; 
residual_R = Rin*Rin' ;
%%
% Rin = Rin - gamma * g_reg ;
Rin = Rin - gamma * g ;
%%
logdata(ii).f           = f;
logdata(ii).f_reg       = f_reg;
logdata(ii).g           = g;
logdata(ii).g_reg       = g_reg;
logdata(ii).residual_R  = residual_R;
logdata(ii).Rin         = Rin;
logdata(ii).R_opt       = R - Rin;
logdata(ii).dir_diff    = sign(R-Rin) - sign(g);
%%
end
%%
E = ones(1,6);
f = slra_mex_obj('mhess', obj, Rin, E)
%%
E = [0.1 0.2 0.5 1 2 3];
f = slra_mex_obj('mhess', obj, Rin, E)
%%
f = slra_mex_obj('func', obj, Rin)
g = slra_mex_obj('grad', obj, Rin)
mhess = slra_mex_obj('mhess', obj, Rin, E)

%% DELETE OBJ
slra_mex_obj('delete', obj);
