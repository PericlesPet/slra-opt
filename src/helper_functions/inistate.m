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
