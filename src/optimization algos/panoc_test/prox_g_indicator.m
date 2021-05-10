function [ prox_g_x ] = prox_g_indicator( R )

% If R is already stiefel, dont waste time projecting again
    if g_indicator(R) <= 1
        prox_g_x = R(:);
% Else, project onto the stiefel manifold
    else
        projR = stiefSVD_proj(R);
        prox_g_x = projR(:);
    end
    
end

