function th = R2th(R, d, psi, R0)
if (size(psi, 1) ~= size(psi, 2)) && (all(all(R0 == 0)))
  P = null(R); dh = P * (P \ d); size(dh);
  th = lra(psi * kron(dh, eye(size(R, 1))), size(psi, 1) - 1);
else
  th = vec(R - R0)' / psi;
end
