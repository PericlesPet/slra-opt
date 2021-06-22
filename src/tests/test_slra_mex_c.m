% USAGE
%     tic
%     [ph_test, info_test, wtfdata2] = test_slra_mex_c(w, m_in, ell, opt_oe);
%     t_test = toc;
%     fprintf('test_slra_mex_c took %f seconds', t_test);
%     
%     wtfdata2.w_ini = w;
%     wtfdata2.opt_ini = opt_oe;
%     wtfdata2.p_ini = w2p(w);
%     wtfdata2.s_ini = s;
function [ph, info, wtfdata2] = test_slra_mex_c(w, m, ell, opt)

wtfdata2.w      = w;
wtfdata2.opt    = opt;

[w, s, r, opt, q, N, T] = sys2slra(w, m, ell, opt);
p = w2p(w);

wtfdata2.w_new = w;
wtfdata2.s_new = s;
wtfdata2.opt_new = opt;
wtfdata2.p      = p;



for opt_options_handling = 1    
if ~exist('opt'), opt = struct; end
if ~isfield(opt, 'solver'), opt.solver = 'c'; end 
if ~isfield(opt, 'disp'), opt.disp = 'off'; end 
if ~isfield(opt, 'tol_m'), opt.tol_m = 1e-6; end 
end 

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

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% IF SOLVER = C %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.solver == 'c'
  if isfield(s, 'w'), s.w(find(s.w == 0)) = opt.tol_m; end
  opt = rmfield(opt, 'solver');
  obj = slra_mex_obj('new', p, s, r);
  [ph, info] = slra_mex_obj('optimize', obj, opt);
  slra_mex_obj('delete', obj);
end







function S = s2s(s, np)
q = length(s.m); if exist('p'), 
                   np = length(p); 
                 else
                   np = sum(s.m) * length(s.n) + length(s.m) * sum(s.n) ...
                                               - length(s.m) * length(s.n);
                 end, if ~isfield(s, 'n'), s.n = np - sum(s.m) + 1; end
N = length(s.n); n  = sum(s.n);, p = 1:np;  
if ~isfield(s, 'phi'), s.phi = eye(sum(s.m)); end, [m, mp] = size(s.phi); 
tmp = cumsum([1; s.m(:)]); Imb = tmp(1:end - 1); Ime = tmp(2:end) - 1;
tmp = cumsum([1; s.n(:)]); Inb = tmp(1:end - 1); Ine = tmp(2:end) - 1;
S = zeros(mp, n); ind = 1;
for j = 1:N
  for i = 1:q
    npij = s.m(i) + s.n(j) - 1;
    pij = p(ind:(ind + npij - 1)); ind = ind + npij;
    Hij = hankel(pij(1:s.m(i)), pij(s.m(i):end));
    S(Imb(i):Ime(i), Inb(j):Ine(j)) = Hij;
  end
end

