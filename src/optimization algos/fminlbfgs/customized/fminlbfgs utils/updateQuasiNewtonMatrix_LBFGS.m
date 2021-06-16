function data = updateQuasiNewtonMatrix_LBFGS(data,optim)
% updates the quasi-Newton matrix that approximates the inverse to the Hessian.
% Two methods are support BFGS and L-BFGS, in L-BFGS the hessian is not
% constructed or stored.
% Calculate position, and gradient diference between the
% itterations
deltaX=data.alpha* data.dir;
deltaG=data.gradient-data.gOld;
        
if ((deltaX'*deltaG) >= sqrt(eps)*max( eps,norm(deltaX)*norm(deltaG) ))
    if(optim.HessUpdate(1)=='b')
        % Default BFGS as described by Nocedal
        p_k = 1 / (deltaG'*deltaX);
        Vk = eye(data.numberOfVariables) - p_k*deltaG*deltaX';
        % Set Hessian
        data.Hessian = Vk'*data.Hessian *Vk + p_k * deltaX*deltaX'
        % Set new Direction
        data.dir = -data.Hessian*data.gradient;
    else
        % L-BFGS with scaling as described by Nocedal
       
        % Update a list with the history of deltaX and deltaG
        data.deltaX(:,2:optim.StoreN)=data.deltaX(:,1:optim.StoreN-1); data.deltaX(:,1)=deltaX;
        data.deltaG(:,2:optim.StoreN)=data.deltaG(:,1:optim.StoreN-1); data.deltaG(:,1)=deltaG;
    
        data.nStored=data.nStored+1; if(data.nStored>optim.StoreN), data.nStored=optim.StoreN; end
        % Initialize variables
        a=zeros(1,data.nStored);
        p=zeros(1,data.nStored);
        q = data.gradient;
        for i=1:data.nStored
            p(i)= 1 / (data.deltaG(:,i)'*data.deltaX(:,i));
            a(i) = p(i)* data.deltaX(:,i)' * q;
            q = q - a(i) * data.deltaG(:,i);
        end
        % Scaling of initial Hessian (identity matrix)
        p_k = data.deltaG(:,1)'*data.deltaX(:,1) / sum(data.deltaG(:,1).^2); 
        
        % Make r = - Hessian * gradient
        r = p_k * q;
        for i=data.nStored:-1:1,
            b = p(i) * data.deltaG(:,i)' * r;
            r = r + data.deltaX(:,i)*(a(i)-b);
        end
        
        % Set new direction
        data.dir = -r;
    end
end
