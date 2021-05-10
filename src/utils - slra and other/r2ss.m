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
