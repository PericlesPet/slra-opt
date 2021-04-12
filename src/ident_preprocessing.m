function [w, s, r, opt, q, N, T] = ident_preprocessing(w,m,ell,opt)

    if ~exist('opt'), opt = []; end
    if ~isfield(opt, 'exct'), opt.exct = []; end
    if ~iscell(w)
      [T, q, N] = size(w); T = ones(N, 1) * T;
    else
      N = length(w); for k = 1:N, [T(k), q, NN(k)] = size(w{k}); end, T = T(:);
      if any(NN > 1), error('W can not be a cell array with a 3D element.'), end
    end, p = q - m; n = p * ell; 

    % IF M == 0
    for M_Equals_Zero = 1
    if (m == 0) & (p > 1)
      if ~iscell(w)
        for k = 1:N, for i = 1:p, wp{(k - 1) * p + i}  = w(:, i, k); end, end
      else
        for k = 1:N, for i = 1:p, wp{(k - 1) * p + i}  = w{k}(:, i); end, end
      end

      if isfield(opt, 'n'), 
        ellp = opt.n; opt = rmfield(opt, 'n'); 
      else, 
        ellp = ell * p; 
      end 

      optp = opt;
      % exact variables TODO
      % initial conditions TODO
      % initial system
      if isfield(optp, 'sys0') % && ???
        optp.sys0 = ss(optp.sys0.a, [], ones(1, size(optp.sys0, 'order')), [], -1);
      end

      [sysph, infop, wph] = ident_custom(wp, 0, ellp, optp); 
      for i = 1:p, c(p, :) = inistate(wph{i}, sysph)'; end
      sysh = ss(sysph.a, [], c, [], -1); 

      if ~iscell(w)
        for k = 1:N, for i = 1:p, wh(:, i, k) = wph{(k - 1) * p + i}; end, end  
      else
        for k = 1:N, for i = 1:p, wh{k}(:, i) = wph{(k - 1) * p + i}; end, end      
      end

      info = infop; % is this correct?
      return 
    end
    end

    % STRUCTURE M, N, r
    s.m = (ell + 1) * ones(q, 1); 
    s.n = T - ell; 
    r = (ell + 1) * m + n;
    %%%%%%%%%%%%%%%%%%%%%%%
    %%% IF EXACT VALUES %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(opt.exct)
      if iscell(opt.exct)
        for i = 1:N, s.w(:, i) = ones(q, 1); s.w(opt.exct{i}) = inf; end
      else 
        s.w = ones(q, 1); s.w(opt.exct) = inf; 
      end
    end

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% IF W_ini in Opt %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    for w_initial_exists = 1
    if isfield(opt, 'wini')
      if iscell(w)
        if isempty(opt.wini) 
          opt.wini = cell(1, N); 
        elseif ~iscell(opt.wini) & (opt.wini == 0)
          for k = 1:N, wini{k} = zeros(ell, q); end, opt.wini = wini;
        else  
          for k = 1:N 
            if opt.wini{k} == 0, opt.wini{k} = zeros(ell, q); end
          end 
        end
      elseif opt.wini == 0
        opt.wini = zeros(ell, q, N);
      end
      if isfield(s, 'w'), s = rmfield(s, 'w'); end
      if ~iscell(opt.wini) && ~isempty(opt.wini)
        W = ones(T, q, N); W(:, opt.exct, :) = inf;
        s.n = s.n + ell; T = T + ell;
        s.w = [inf * ones(ell, q, N); W]; 
        w   = [opt.wini; w];
      elseif iscell(opt.wini) 
        for k = 1:N
          W = ones(T(k), q); W(:, opt.exct{k}) = inf;
          if ~isempty(opt.wini{k})
            s.n(k) = s.n(k) + ell; T(k) = T(k) + ell;
            s.w{k} = [inf * ones(ell, q); W]; 
            w{k}   = [opt.wini{k}; w{k}];
          else
            s.w{k} = W;
          end
        end
      end
    end
    end

    for sys_initial_exists = 1
    if isfield(opt, 'sys0')
      if isa(opt.sys0, 'lti'),  
        opt.Rini = ss2r(opt.sys0); 
      else, 
        opt.Rini = opt.sys0; 
      end
    end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% IF S (structure) HAS FIELD W (weights) %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(s, 'w'),
      if isfield(opt, 'wini'), s.w = w2p(s.w); end
      if all(size(s.w) == [q 1]), s.w = s.w(:, ones(1, N)); end
    end
    if isfield(opt, 'v'),
      if all(size(opt.v) == [q 1]), opt.v = opt.v(:, ones(1, N)); end
      if ~all(size(opt.v) == [q N]), opt.v = w2p(opt.v); end
    end
    if isfield(s, 'w') && isfield(opt, 'v')
      try, s.w = s.w .* opt.v; catch
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
        w_inf = s.w; s.w = opt.v; 
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
        s.w = s.w .* w_inf;
      end
    elseif isfield(opt, 'v'), s.w = opt.v; end
