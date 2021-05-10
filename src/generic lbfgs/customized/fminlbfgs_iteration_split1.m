function [x,fval,exitflag,output,grad, data]=fminlbfgs_iteration_split1(funfcn, optim, data, exitflag, hessian_update)
%   Function is written by D.Kroon University of Twente (Updated Nov. 2010)
% Removed comments -> add them in final version     

% Start Minimizing
while(true)
    
    if (~hessian_update) 
    % Update number of itterations
    data.iteration=data.iteration+1; 
    % Set current lineSearch parameters
    data.TolFunLnS = eps(max(1,abs(data.fInitial )));
    data.fminimum = data.fInitial - 1e16*(1+abs(data.fInitial));
    
	% Make arrays to store linesearch results
    data.storefx=[]; data.storepx=[]; data.storex=[]; data.storegx=[];
    % If option display plot, than start new figure
    if(optim.Display(1)=='p'), figure, hold on; end
		
    % Find a good step size in the direction of the gradient: Linesearch
    if(optim.GoalsExactAchieve==1)
		data=linesearch(funfcn, data,optim);
    else
        data=linesearch_simple(funfcn, data, optim);
    end
	
	% Make linesearch plot
	if(optim.Display(1)=='p'); 
		plot(data.storex,data.storefx,'r*');
		plot(data.storex,data.storefx,'b');
		
		alpha_test= linspace(min(data.storex(:))/3, max(data.storex(:))*1.3, 10);
		falpha_test=zeros(1,length(alpha_test));
        for i=1:length(alpha_test)
			[data,falpha_test(i)]=gradient_function(data.xInitial(:)+alpha_test(i)*data.dir(:),funfcn, data, optim);
        end    
		plot(alpha_test,falpha_test,'g');
        plot(data.alpha,data.f_alpha,'go','MarkerSize',8);
	end
	
    % Check if exitflag is set
    if(~isempty(data.exitflag)),
        exitflag=data.exitflag;
        data.xInitial=data.xOld; 
        data.fInitial=data.fOld;
        data.gradient=data.gOld;
        break, 
    end;
    
    % Update x with the alpha step
    data.xInitial = data.xInitial + data.alpha*data.dir;
    
    % Set the current error and gradient
    data.fInitial =  data.f_alpha;
	data.gradient = data.grad;
    
    % Set initial steplength to 1
    data.initialStepLength = 1;
    
    
    gNorm = norm(data.gradient,Inf);  % Norm of gradient
    
    % Set exit flags 
    if(gNorm <optim.TolFun), exitflag=1; end
    if(max(abs(data.xOld-data.xInitial)) <optim.TolX), exitflag=2; end
    if(data.iteration>=optim.MaxIter), exitflag=0; end
    
    % Check if exitflag is set
    if(~isempty(exitflag)), break, end;
    
%   STOP BEFORE HESSIAN UPDATE
    else
        
    % Update the inverse Hessian matrix
    if(optim.HessUpdate(1)~='s')
        % Do the Quasi-Neton Hessian update.
        data = updateQuasiNewtonMatrix_LBFGS(data,optim);
    else
        data.dir = -data.gradient;
    end
  
    % Derivative of direction
    data.fPrimeInitial= data.gradient'*data.dir(:);
    % Call output function
    if(call_output_function(data,optim,'iter')), exitflag=-1; end
    
    % Show the current iteration
    if(strcmp(optim.Display(1),'i')||strcmp(optim.Display(1),'p'))
        s=sprintf('     %5.0f       %5.0f       %5.0f       %13.6g   %13.6g',data.iteration,data.funcCount,data.gradCount,data.fInitial,data.alpha); disp(s);
    end
    
    % Keep the variables for next iteration
    data.fOld=data.fInitial;
    data.xOld=data.xInitial;
    data.gOld=data.gradient;
    end
    break;
end
% Set output parameters
fval=data.fInitial;
grad=data.gradient;
x = data.xInitial;
% Reshape x to original shape
x=reshape(x,data.xsizes);
% Call output function
if(call_output_function(data,optim,'done')), exitflag=-1; end

% Make exist output structure
if ~isempty(exitflag)
if(optim.HessUpdate(1)=='b'), output.algorithm='Broyden?letcher?oldfarb?hanno (BFGS)';
elseif(optim.HessUpdate(1)=='l'), output.algorithm='limited memory BFGS (L-BFGS)';
else output.algorithm='Steepest Gradient Descent'; 
end
output.message=getexitmessage(exitflag);
output.iteration = data.iteration;
output.funccount = data.funcCount;
output.fval = data.fInitial;
output.stepsize = data.alpha;
output.directionalderivative = data.fPrimeInitial;
output.gradient = reshape(data.gradient, data.xsizes);
output.searchdirection = data.dir;
output.timeTotal=toc(data.timeTotal);    
output.timeExtern=data.timeExtern;
oupput.timeIntern=output.timeTotal-output.timeExtern;
% Display final results
if(~strcmp(optim.Display,'off'))
    disp('    Optimizer Results')
    disp(['        Algorithm Used: ' output.algorithm]);
    disp(['        Exit message : ' output.message]);
    disp(['        iterations : '  int2str(data.iteration)]);
    disp(['        Function Count : ' int2str(data.funcCount)]);
    disp(['        Minimum found : ' num2str(fval)]);
    disp(['        Intern Time : ' num2str(oupput.timeIntern) ' seconds']);
    disp(['        Total Time : ' num2str(output.timeTotal) ' seconds']);
end
else
    output = null(1);
end