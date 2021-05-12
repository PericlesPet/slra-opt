% THIS FUNCTION RETURNS THE VALUES OF F(X) AND dF(X) 
% DOES NOT RETURN HANDLES
% CAN BE USED FOR ANONYMOUS DECLARATIONS 
% E.G. problem_func = @(x) SLRA_FuncAndGrad_vals(obj, mu, dimensions, x, mode)
%      [f_val, df_val] = problem_func(x);
function [fVal, dfVal] = SLRA_FuncAndGrad_vals(obj, mu, dimensions, x, fcn_mode)

if nargin<5, fcn_mode = []; end
[f, df] = SLRA_FuncAndGrad_handles(obj, mu, dimensions, fcn_mode);


fVal    = f(x);
if nargout > 1
    dfVal   = df(x);    
end

end

% OLD CODE

% dimension = dimensions(1) * dimensions(2);
% 
% % Rkern is the Kernel of the Hankel matrix
% Rkern = @(x) reshape(x, dimensions);
% 
% 
% % ASSIGN F
% f   = @(x) slra_mex_obj('func', obj, Rkern(x));
% if nargin == 5     % IF A MODE IS SELECTED
%     if strcmp(mode, 'reg')
%         f =@(x) f(x) + mu * norm(Rkern(x) * Rkern(x)' - ...
%              eye(dimensions(1)),'fro')^2;
%     elseif strcmp(mode, 'alm')    
%         %%TO DO
%     end
% end
% fVal    = f(x);
% 
% % ASSIGN DF IF NECESSARY
% if nargout > 1
%     df  = @(x) reshape(slra_mex_obj('grad', obj, Rkern(x)), dimension , 1);
%     
%     if nargin == 5     % IF A MODE IS SELECTED
%         if strcmp(mode, 'reg')
%             grad_regularizer = @(x)  reshape(2 *(Rkern(x)* Rkern(x)' ... 
%                 - eye(dimensions(1)))*Rkern(x), ... 
%                 dimension, 1);
%             df = @(x) df(x) + mu * grad_regularizer(x) ; 
%         elseif strcmp(mode, 'alm')    
%             %%TO DO
%         end            
%     end
%     dfVal   = df(x);    
% end


