
for initial_stuff_of_slraext = 1
    
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


prob.options = optimset(opt); 
if isempty(prob.options.Display)
   prob.options = optimset(prob.options, 'disp', 'off'); 
end

global Th F; Th = []; F = []; 
prob.options = optimset(prob.options, 'OutputFcn', @OutputFcn);

prob.solver = 'fmincon'; 

pext = [0; p]; 

prob.x0 = R2th(Rini, phi * (s0 + pext(tts + 1)), psi, opt.R0); 

%%
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


%%
prob.objective = @(th) slra_mex_obj('func', obj, th2R(th));

prob.nonlcon = @(th) deal([], C(th));
prob.options = optimset(prob.options, 'alg', 'sqp');
[x, fval, flag, info] = fmincon(prob); 
%%

info.fmin = fval;, info.Rh = th2R(x); 
info.Th = Th; info.F = F;
% [M, ph] = Mslra_ext(info.Rh, tts, p, [], bfs, phi, s0);
[M, ph] = Mslra_ext(info.Rh, tts, p, [], bfs, phi, s0);

M2 = slra_mex_obj('func', obj, info.Rh)
%%
if exist('w') & ~isempty(w) 
  ph(Ig) = inv_sqrt_w * ph(Ig); 
  if ~isempty(If)
    ph_ = Inf * ones(np_, 1);
    ph_(Ig_(If)) = pf; ph_(isinf(ph_)) = ph;
    ph = ph_;
  end  
end