function x = almSearch(L, dL, x0, lambda, rho, omega, lb, ub, Xtol, alpha, miter, slradata)
    % ALMSEARCH Find a local minimum of a function using projected gradient descent and
    %           and golden section line search.
    %
    %   Dependencies: almProj.m
    %
    %   Written by Joel T. Kaardal, July 7th, 2016

    MAXITS = miter;
    GOLD = (1+sqrt(5))/2;
    ALPHA = alpha;
    TOL = max([Xtol, omega]);
%     TOL = 0.2;
    grad = dL(x0, lambda, rho);
    x = almProj(x0, ALPHA*grad, lb, ub);
    its = 1;
    iteration_num = 1;
%     while norm(x-x0) > TOL && iteration_num <= 12
    while norm(x-x0) > 0.1
        iteration_num = iteration_num + 1;
        % projected gradient descent with projected line search
        a = x0;
        c = almProj(x0, ALPHA*grad*(GOLD-1)/GOLD, lb, ub);
        d = almProj(x0, ALPHA*grad/GOLD, lb, ub);
        GAMMA = ALPHA;
        while norm(c-d) > TOL/2
            if L(c, lambda, rho) > L(d, lambda, rho)
                a = c;
            end
            GAMMA = GAMMA/GOLD;

            grad = dL(a, lambda, rho);
            c = almProj(a, GAMMA*grad*(GOLD-1)/GOLD, lb, ub);
            d = almProj(a, GAMMA*grad/GOLD, lb, ub);
        end

        x = (a+almProj(a, GAMMA*grad, lb, ub))/2;
        if its > MAXITS
            fprintf('Maximum iterations exceeded in almSearch.\n');
            break;
        end
        its = its+1;
        if norm(x-x0) > TOL
            x0 = x;
            grad = dL(x0, lambda, rho);
            %BETA = grad'*grad/(grad_old'*grad_old);
            x = almProj(x0, ALPHA*grad, lb, ub);
        end
        if mod(iteration_num, 1) == 0
            R_alm0  = reshape(x0(slradata.np+1:end), size(slradata.Rini));
            R_alm   = reshape(x(slradata.np+1:end), size(slradata.Rini));

            fprintf('   INNER ITERATION = %d, norm(x-x0) = %f\n', ... 
                iteration_num, ...
                norm(x-x0));
            fprintf('   	f = %4.4f,         f0 = %4.4f\n', ... 
                slra_mex_obj('func', slradata.obj, R_alm), ...
                slra_mex_obj('func', slradata.obj, R_alm0));
            fprintf('       L = %4.4f,         L0 = %4.4f\n', ... 
                L(x, lambda, rho), ...
                L(x0, lambda, rho));
            fprintf('       CE = %4.4f,        CE0 = %4.4f\n', ... 
                slradata.ce(x), ...
                slradata.ce(x0));
            fprintf('       DP = %4.4f,        DP0 = %4.4f\n', ... 
                norm(x(1:slradata.np) - slradata.p), ... 
                norm(x0(1:slradata.np) - slradata.p));
%             u0 = 
%             sys_comparison(u0(1:length(u0)),:),y0(1:ceil(length(y0)/sample_divisor1),:), r2ss(R_ident, m_in, ell));
%             [Yh, M] = compare(iddata(y, u), idss(sys_ID)); 
        end
    end
end
