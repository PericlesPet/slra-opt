function C = stiefConstraint(R, metric)

    if nargin == 1  % Default: R * R'
        C = R * R';

    elseif nargin > 1 && strcmp(metric,'dist')
    % 'Distance' of R from Stiefel Manifold
        C = norm(R*R' - eye(size(R,1)));

    elseif nargin > 1 && strcmp(metric,'matrix')
    % Otherwise return the Constraint Matrix    
        C = R*R' - eye(size(R,1));
    end

end