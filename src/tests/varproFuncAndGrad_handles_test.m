function [f, df] = varproFuncAndGrad_handles_test(obj, mu, dimensions, mode)

dimension = dimensions(1) * dimensions(2);

% Rkern is the Kernel of the Hankel matrix
Rkern = @(x) reshape(x, dimensions);


% ASSIGN F
f   = @(x) slra_mex_obj('func', obj, Rkern(x));
if nargin == 4     % IF A MODE IS SELECTED
    if strcmp(mode, 'reg')
        f =@(x) f(x) + mu * norm(Rkern(x) * Rkern(x)' - ...
             eye(dimensions(1)),'fro')^2;
    elseif strcmp(mode, 'alm')    
        %%TO DO
    end
end

% ASSIGN DF IF NECESSARY
if nargout > 1
    df  = @(x) reshape(slra_mex_obj('grad', obj, Rkern(x)), dimension , 1);
    
    if nargin == 4     % IF A MODE IS SELECTED
        if strcmp(mode, 'reg')
            grad_regularizer = @(x)  reshape(2 *(Rkern(x)* Rkern(x)' ... 
                - eye(dimensions(1)))*Rkern(x), ... 
                dimension, 1);
            df = @(x) df(x) + mu * grad_regularizer(x) ; 
        elseif strcmp(mode, 'alm')    
            %%TO DO
        end            
    end
end

end

