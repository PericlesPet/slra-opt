% Transform
% FROM: slra(p, s, r, opt) 
% TO:   slra_ext(tts, p, r, w, Rini, phi, psi, opt, th2R, C, s0) ARGUMENTS
function [tts, p, r, w, Rini, phi ,psi, opt, th2R, C, s0] = slra_to_slraext(p, s, r, opt)

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

for if_solver_equals_M = 1
if ~isfield(s, 'w'), s.w = []; end
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
if ~all(size(s.w) == [length(p(:)) length(p(:))]), s.w = s.w(:); end
if ~isfield(s,   'phi' ), s.phi = []; end
if ~isfield(opt, 'Rini'), opt.Rini = []; end
if ~isfield(opt, 'psi' ), opt.psi = []; end
np = length(p); warning_state = warning; warning('off');


tts  = s2s(s, np);
p    = p ;
r    = r;
w    = s.w ;
Rini = opt.Rini;
phi  = s.phi ;
psi  = opt.psi;
opt  = opt;

end
  
for if_exists_in_slra_ext = 1 

p = p(:); if ~exist('opt'), opt = []; end
[mp, n] = size(tts); np = max(max(tts));
if exist('phi', 'var') && ~isempty(phi), m = size(phi, 1); else m = mp; end 
vec_tts = tts(:); NP = 1:np;
bfs = vec_tts(:, ones(1, np)) == NP(ones(mp * n, 1), :);
if ~exist('phi', 'var') | isempty(phi), phi = eye(size(tts, 1)); end
if ~exist('s0') || isempty(s0), s0 = zeros(mp, n); end
if ~exist('psi', 'var') | isempty(psi), psi = eye(m * (m - r)); end
if ~isfield(opt, 'R0'), opt.R0 = 0; end
if ~exist('th2R') || isempty(th2R) 
  th2R = @(th) reshape(th * psi + opt.R0(:)', m - r, m); 
end 
if ~exist('C') || isempty(C), C = @(th) th2R(th) * th2R(th)' - eye(m - r); end  
if ~exist('Rini') | isempty(Rini)
  pext = [0; p];
  Rini = lra(phi * (s0 + pext(tts + 1)), r);
end
prob.options = optimset(opt); 
if isempty(prob.options.Display)
   prob.options = optimset(prob.options, 'disp', 'off'); 
end
    
end