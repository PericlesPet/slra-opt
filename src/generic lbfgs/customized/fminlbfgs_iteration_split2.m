function [x,fval,exitflag,output,grad, data]=fminlbfgs_iteration_split2(funfcn, optim, data, exitflag)
%FMINLBFGS finds a local minimum of a function of several variables. 
%   This optimizer is developed for image registration methods with large 
%	amounts of unknown variables.
%
%   Optimization methods supported:
%	- Quasi Newton Broyden?letcher?oldfarb?hanno (BFGS)  
%   - Limited memory BFGS (L-BFGS)
%   - Steepest Gradient Descent optimization.
%   
%   [X,FVAL,EXITFLAG,OUTPUT,GRAD] = FMINLBFGS(FUN,X0,OPTIONS) 
%
%   Inputs,
%		FUN: Function handle or string which is minimized, returning an
%				error value and optional the error gradient. 
%		X0: Initial values of unknowns can be a scalar, vector or matrix
%	 (optional)
%		OPTIONS: Structure with optimizer options, made by a struct or
%				optimset. (optimset doesnot support all input options)
%
%   Outputs,
%		X : The found location (values) which minimize the function.
%		FVAL : The minimum found
%		EXITFLAG : Gives value, which explain why the minimizer stopt
%		OUTPUT : Structure with all important ouput values and parameters
%		GRAD : The gradient at this location 
%
%   Extended description of input/ouput variables 
%   OPTIONS,
%		OPTIONS.GoalsExactAchieve : If set to 0, a line search method is
%               used which uses a few function calls to do a good line
%               search. When set to 1 a normal line search method with Wolfe 
%				conditions is used (default).
%		OPTIONS.GradConstr, Set this variable to true if gradient calls are
%				cpu-expensive (default). If false more gradient calls are 
%				used and less function calls.
%	    OPTIONS.HessUpdate : If set to 'bfgs', Broyden?letcher?oldfarb?hanno 
%				optimization is used (default), when the number of unknowns is 
%				larger then 3000 the function will switch to Limited memory BFGS, 
%				or if you set it to 'lbfgs'. When set to 'steepdesc', steepest 
%				decent optimization is used.
%		OPTIONS.StoreN : Number of itterations used to approximate the Hessian,
%			 	in L-BFGS, 20 is default. A lower value may work better with
%				non smooth functions, because than the Hessian is only valid for
%				a specific position. A higher value is recommend with quadratic equations. 
%		OPTIONS.GradObj : Set to 'on' if gradient available otherwise finited difference
%				is used.
%     	OPTIONS.Display : Level of display. 'off' displays no output; 'plot' displays
%				all linesearch results in figures. 'iter' displays output at  each 
%               iteration; 'final' displays just the final output; 'notify' 
%				displays output only if the function does not converge; 
%	    OPTIONS.TolX : Termination tolerance on x, default 1e-6.
%	    OPTIONS.TolFun : Termination tolerance on the function value, default 1e-6.
%		OPTIONS.MaxIter : Maximum number of iterations allowed, default 400.
% 		OPTIONS.MaxFunEvals : Maximum number of function evaluations allowed, 
%				default 100 times the amount of unknowns.
%		OPTIONS.DiffMaxChange : Maximum stepsize used for finite difference gradients.
%		OPTIONS.DiffMinChange : Minimum stepsize used for finite difference gradients.
%		OPTIONS.OutputFcn : User-defined function that an optimization function calls
%				at each iteration.
%		OPTIONS.rho : Wolfe condition on gradient (c1 on wikipedia), default 0.01.
%		OPTIONS.sigma : Wolfe condition on gradient (c2 on wikipedia), default 0.9. 
%		OPTIONS.tau1 : Bracket expansion if stepsize becomes larger, default 3.
%		OPTIONS.tau2 : Left bracket reduction used in section phase,
%		default 0.1.
%		OPTIONS.tau3 : Right bracket reduction used in section phase, default 0.5.
%   FUN,
%		The speed of this optimizer can be improved by also providing
%   	the gradient at X. Write the FUN function as follows
%   	function [f,g]=FUN(X)
%       	f , value calculation at X;
%   	if ( nargout > 1 )
%       	g , gradient calculation at X;
%   	end
%	EXITFLAG,
%		Possible values of EXITFLAG, and the corresponding exit conditions
%		are
%  		1, 'Change in the objective function value was less than the specified tolerance TolFun.';
%  		2, 'Change in x was smaller than the specified tolerance TolX.'; 
%  		3, 'Magnitude of gradient smaller than the specified tolerance';
%  		4, 'Boundary fminimum reached.';
%  		0, 'Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
%  		-1, 'Algorithm was terminated by the output function.';
%  		-2, 'Line search cannot find an acceptable point along the current search';
%
%   Examples
%       options = optimset('GradObj','on');
%       X = fminlbfgs(@myfun,2,options)
%
%   	% where myfun is a MATLAB function such as:
%       function [f,g] = myfun(x)
%       f = sin(x) + 3;
%	    if ( nargout > 1 ), g = cos(x); end
%
%   See also OPTIMSET, FMINSEARCH, FMINBND, FMINCON, FMINUNC, @, INLINE.
%
%   Function is written by D.Kroon University of Twente (Updated Nov. 2010)

% Start Minimizing
while(true)
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





    
    
    