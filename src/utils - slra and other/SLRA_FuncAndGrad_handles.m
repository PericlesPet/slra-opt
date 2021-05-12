function [f, df] = SLRA_FuncAndGrad_handles(obj, mu, dimensions, fcn_mode)

dimension = dimensions(1) * dimensions(2);

% Rkern is the Kernel of the Hankel matrix
Rkern = @(x) reshape(x, dimensions);


% ASSIGN F
f   = @(x) slra_mex_obj('func', obj, Rkern(x));
if nargin == 4     % IF A MODE IS SELECTED
    if strcmp(fcn_mode, 'reg')
        f =@(x) f(x) + mu * norm(Rkern(x) * Rkern(x)' - ...
             eye(dimensions(1)),'fro')^2;
    elseif strcmp(fcn_mode, 'alm')    
        %%TO DO
    end
end

% ASSIGN DF IF NECESSARY
if nargout > 1
    df  = @(x) reshape(slra_mex_obj('grad', obj, Rkern(x)), dimension , 1);
    
    if nargin == 4     % IF A MODE IS SELECTED
        if strcmp(fcn_mode, 'reg')
            grad_regularizer = @(x)  reshape(2 *(Rkern(x)* Rkern(x)' ... 
                - eye(dimensions(1)))*Rkern(x), ... 
                dimension, 1);
            df = @(x) df(x) + mu * grad_regularizer(x) ; 
        elseif strcmp(fcn_mode, 'alm')    
            %%TO DO
        end            
    end
%   SCALING
%     scale_factor = 1/1000;
    alpha = 0.5;
    scale_factor = @(x) alpha * norm(x) / norm(df(x));
    df = @(x) scale_factor(x) * df(x);
end

end

