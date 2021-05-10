function [data,fval,grad]=gradient_function(x,funfcn, data, optim)
    % Call the error function for error (and gradient)
    if ( nargout <3 )
        timem=tic;   
        fval=funfcn(reshape(x,data.xsizes)); 
        data.timeExtern=data.timeExtern+toc(timem);
        data.funcCount=data.funcCount+1;
    else
        if(strcmp(optim.GradObj,'on'))
            timem=tic;    
            [fval, grad]=feval(funfcn,reshape(x,data.xsizes)); 
            data.timeExtern=data.timeExtern+toc(timem);
            data.funcCount=data.funcCount+1;
            data.gradCount=data.gradCount+1;
        else
            % Calculate gradient with forward difference if not provided by the function
            grad=zeros(length(x),1);
            fval=funfcn(reshape(x,data.xsizes));
            gstep=data.initialStepLength/1e6; 
            if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
            if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
            for i=1:length(x),
                x_temp=x; x_temp(i)=x_temp(i)+gstep;
                timem=tic;    
                [fval_g]=feval(funfcn,reshape(x_temp,data.xsizes)); data.funcCount=data.funcCount+1;
                data.timeExtern=data.timeExtern+toc(timem);
                grad(i)=(fval_g-fval)/gstep;
            end
        end
        grad=grad(:);
    end
