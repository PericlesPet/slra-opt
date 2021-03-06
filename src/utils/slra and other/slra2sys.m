function [sysh, info, wh, xini] = slra2sys(w, m, ell, opt, ph, info, q, N, T)

wh = p2w(ph, q, N, T, iscell(w)); 
if isfield(opt, 'ss') && (opt.ss == 0) 
  sysh = info.Rh; xini = NaN;
else 
  sysh = r2ss(info.Rh, m, ell); 
  xini = inistate(wh, sysh);
end
if isfield(opt, 'wini')
  if ~iscell(opt.wini) && ~isempty(opt.wini)
    wh = wh(ell + 1:end, :, :);
  elseif iscell(opt.wini) && ~isempty(cell2mat(opt.wini))
    for k = 1:N
      if ~isempty(opt.wini{k})
        wh{k} = wh{k}(ell + 1:end, :);
      end
    end  
  end
end
if isfield(opt, 'n') && (p * ell) ~= opt.n
  sysh = balred(sysh, opt.n);
end

% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   W2P   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = w2p(w)
if ~iscell(w), p = w(:); else
  p = []; 
  for k = 1:length(w), p = [p; w2p(w{k})]; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   P2W   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function w = p2w(p, q, N, T, c)
if ~c 
  if N == 1, w = reshape(p, T, q, N); else
    for k = 1:N, w(:, :, k) = p2w(p(1:q * T(k)), q, 1, T(k), 0); p(1:q * T(k)) = []; end
  end
else
  for k = 1:N, w{k} = p2w(p(1:q * T(k)), q, 1, T(k), 0); p(1:q * T(k)) = []; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   SS2R  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = ss2r(sys)
sys = ss(sys); a = sys.a; b = sys.b; c = sys.c; d = sys.d; 
[p, m] = size(d); n = size(a, 1); 
ell1 = n / p + 1; L = ell1; O = c; for t = 2:L, O = [O; O(end - p + 1:end, :) * a]; end, P = null(O')';
if (m > 0)
  F = [d; O(1:(end - p), :) * b]; TT = zeros(ell1 * p, ell1 * m);
  for i = 1:ell1
    TT((i - 1) * p + 1:end, (i - 1) * m + 1: i * m) = F(1:(ell1 + 1 - i) * p, :);
  end
  Q = P * TT; 
else, Q = []; end
R = permute([reshape(Q, p, m, ell1), -reshape(P, p, p, ell1)], [1 3 2]);
R = reshape(R, p, (m + p) * ell1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   R2SS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sysh = r2ss(R, m, ell)
[p, tmp] = size(R); ell1 = ell + 1; q = tmp / ell1; n = ell * p;
R = permute(reshape(R, p, ell1, q), [1 3 2]);
Q = R(:, 1:m, :); P = - R(:, m + 1:q, :); inv_Pl = pinv(P(:, :, ell + 1));
a = zeros(n); b = zeros(n, m); c = [];
if n > 0
  a(p + 1:end, 1:n - p) = eye(n - p); 
  c = [zeros(p, n - p) eye(p)];
end
d = inv_Pl * Q(:, :, ell1); ind_j = (n - p + 1):n;
for i = 1:ell
  ind_i = ((i - 1) * p + 1):(i * p); Pi = P(:, :, i);
  a(ind_i, ind_j) = - inv_Pl * Pi; 
  b(ind_i, :) = inv_Pl * (Q(:, :, i) - Pi * d);
end
sysh = ss(a, b, c, d, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%   XINI   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xini = inistate(w, sys, use_all_data)
a = sys.a; c = sys.c; [p, m] = size(sys.d); n = size(a, 1); 
if ~iscell(w)
  [T, q, N] = size(w); T = ones(N, 1) * T;
else
  N = length(w); for k = 1:N, [T(k), q, NN(k)] = size(w{k}); end, T = T(:);
  if any(NN > 1), error('W can not be a cell array with a 3D element.'), end
end 
if ~exist('use_all_data') || use_all_data ~= 1, T = max(n, 2) * ones(1, N); end
L = max(T); O = c; for t = 2:L, O = [O; O(end - p + 1:end, :) * a]; end
sys.Ts = -1; xini = zeros(n, N);  
for k = 1:N
  if ~iscell(w)
    uk = w(1:T(k), 1:m, k); yk = w(1:L, (m + 1):end, k);
  else  
    uk = w{k}(1:T(k), 1:m); yk = w{k}(1:L, (m + 1):end);
  end
  if m > 0, y0k = (yk - lsim(sys, uk, 0:(T(k) - 1)))'; else, y0k = yk'; end
  xini(:, k) = O(1:(T(k) * p), :) \ y0k(:);
end
