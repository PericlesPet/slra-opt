function [sysh, info, wh, xini] = ident_custom(w, m, ell, opt)

[w, s, r, opt, q, N, T] = sys2slra(w, m, ell, opt);

tic
%%%%%%%%% SLRA %%%%%%%%%%%%%%
[ph, info] = slra(w2p(w), s, r, opt); info.M = info.fmin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;

[sysh, info, wh, xini] = slra2sys(w, m, ell, opt, ph, info, q, N, T);
