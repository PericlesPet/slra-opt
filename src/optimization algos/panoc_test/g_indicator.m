function [ g_x ] = g_indicator( R )

% C = norm(R'*R - I), or "how much" the constraint is violated
    C = stiefConstraint(R, 'dist');
% The indicator function should be inf or 0
% so for simplicity multiply by 10^10
    g_x = C*10^10;

end

