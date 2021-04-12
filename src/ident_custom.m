function [sysh, info, wh, xini] = ident_custom(w, m, ell, opt)

[w, s, r, opt, q, N, T] = ident_preprocessing(w, m, ell, opt);

%%%%%%%%% SLRA %%%%%%%%%%%%%%
[ph, info] = slra(w2p(w), s, r, opt); info.M = info.fmin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sysh, info, wh, xini] = ident_postprocessing(w, m, ell, opt, ph, info, q, N, T);
