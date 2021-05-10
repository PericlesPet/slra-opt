% Transform
% FROM: slra(p, s, r, opt) 
% TO:   slra_ext(tts, p, r, w, Rini, phi, psi, opt, th2R, C, s0) ARGUMENTS
function [tts, p, r, s, w, Rini, phi ,psi, opt, th2R, C, s0, prob, pext] = slra2slra_ext(p, s, r, opt)

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

for if_solver_equals_x = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% IF SOLVER = C %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt.solver == 'c'
  if isfield(s, 'w'), s.w(find(s.w == 0)) = opt.tol_m; end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% IF SOLVER = R %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif opt.solver == 'r'
  if isfield(s, 'w')
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
    opt.w = s.w; 
  end
  if isfield(opt, 'Rini'), opt.P_init = null(opt.Rini); end
  np = length(p); s.tts = s2s(s, np); 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% IF SOLVER = M %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
else
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
  tts = s2s(s, np);
  w = s.w; 
  Rini = opt.Rini;
  phi = s.phi;
  psi = opt.psi;
    
end
end

  
for if_exists_in_slra_ext = 1 
    p_temporary = p;
if opt.solver == 'm'

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


global Th F; Th = []; F = []; 
prob.options = optimset(prob.options, 'OutputFcn', @OutputFcn);
reg = exist('opt') && isfield(opt, 'method') && strcmp(opt.method, 'reg');
if reg
  prob.solver = 'fminunc';
else
  prob.solver = 'fmincon'; 
end
pext = [0; p]; 
%Rini = - Rini(:, end - size(Rini, 1) + 1:end) \ Rini; % quick and dirty fix needed for slra_armax
prob.x0 = R2th(Rini, phi * (s0 + pext(tts + 1)), psi, opt.R0); 
Im = find(isnan(p)); 
if exist('w') && length(w(:)) == length(p), 
  Im = unique([Im(:); find(w(:) == 0)]); 
  p(Im) = NaN;
end
Ig = setdiff(1:np, Im); 
if exist('w') & ~isempty(w)
  if any(size(w) == 1), w = diag(w); end
  if size(w, 1) == np, w = w(Ig, Ig); end  
  If = find(isinf(diag(w))); 
  if ~isempty(If)
    pf = p(Ig(If));
    s0 = s0 + reshape(bfs(:, Ig(If)) * pf, mp, n);
    w(If, :) = []; w(:, If) = []; p(Ig(If)) = []; 
    bfs(:, Ig(If)) = []; 
    Ig_ = Ig; np_ = np; np = length(p); 
    tts = reshape(bfs * vec(1:np), mp, n);
    Im = find(isnan(p)); 
    if exist('w') && length(w(:)) == length(p), 
      Im = unique([Im(:); find(w(:) == 0)]); 
      p(Im) = NaN;
    end
    Ig = setdiff(1:np, Im); 
  end
  sqrt_w = sqrtm(w); inv_sqrt_w = pinv(sqrt_w);
  bfs = double(bfs); p(Ig) = sqrt_w * p(Ig); 
  bfs(:, Ig) = bfs(:, Ig) * inv_sqrt_w;
end
if reg
  if ~exist('opt') || ~isfield(opt, 'g') || isempty(opt.g) 
    opt.g = norm(p(Ig)) ^ 2; 
  end
  prob.objective = @(th) Mslra_ext(th2R(th), tts, p, [], bfs, phi, s0) ...
                         + opt.g * norm(C(th), 'fro') ^ 2;
else
  prob.objective = @(th) Mslra_ext(th2R(th), tts, p, [], bfs, phi, s0);
%   prob.objective = @(th) Mslra_ext2(th2R(th), tts, p, [], bfs, phi, s0);
  prob.nonlcon = @(th) deal([], C(th));
end  
if reg
%   [x, fval, flag, info] = fminunc(prob);
else
  prob.options = optimset(prob.options, 'alg', 'sqp');
%   [x, fval, flag, info] = fmincon(prob); 
end
% info.fmin = fval;, info.Rh = th2R(x); 
% info.Th = Th; info.F = F;
% [M, ph] = Mslra_ext(info.Rh, tts, p, [], bfs, phi, s0);
% if exist('w') & ~isempty(w) 
%   ph(Ig) = inv_sqrt_w * ph(Ig); 
%   if ~isempty(If)
%     ph_ = Inf * ones(np_, 1);
%     ph_(Ig_(If)) = pf; ph_(isinf(ph_)) = ph;
%     ph = ph_;
%   end  
% end
p = p_temporary;
end



    
end