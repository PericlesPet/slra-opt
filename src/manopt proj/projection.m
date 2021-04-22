% grad = proj(x,grad) 

function Up = projection(X, U)
        
        XtU = multiprod(multitransp(X), U);
        symXtU = multisym(XtU);
        Up = U - multiprod(X, symXtU);
        
% The code above is equivalent to, but faster than, the code below.
%         
%     Up = zeros(size(U));
%     function A = sym(A), A = .5*(A+A'); end
%     for i = 1 : k
%         Xi = X(:, :, i);
%         Ui = U(:, :, i);
%         Up(:, :, i) = Ui - Xi*sym(Xi'*Ui);
%     end

    end