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

