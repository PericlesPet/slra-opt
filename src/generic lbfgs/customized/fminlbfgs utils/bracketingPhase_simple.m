function data = bracketingPhase_simple(funfcn, data,optim)
% Number of itterations
itw=0; 
% Point with smaller value, initial
data.beta=0; 
data.f_beta=data.fInitial; 
data.fPrime_beta=data.fPrimeInitial;
% Initial step is equal to alpha of previous step.
alpha = data.initialStepLength;
% Going up hill
hill=false;
% Search for brackets
while(true)
    % Calculate the error registration gradient
    if(optim.GradConstr)
        [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
        fPrime_alpha=nan;
        grad=nan;
    else
        [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
        fPrime_alpha = grad'*data.dir(:);
    end
    
	% Store values linesearch
	data.storefx=[data.storefx f_alpha]; 
    data.storepx=[data.storepx fPrime_alpha]; 
	data.storex=[data.storex alpha]; 
	data.storegx=[data.storegx grad(:)];
    
    % Update step value
    if(data.f_beta<f_alpha), 
        % Go to smaller stepsize
        alpha=alpha*optim.tau3;
        
        % Set hill variable
        hill=true;
    else
        % Save current minium point
        data.beta=alpha; data.f_beta=f_alpha; data.fPrime_beta=fPrime_alpha; data.grad=grad;
        if(~hill)
            alpha=alpha*optim.tau1;  
        end
    end
                        
    % Update number of loop iterations
    itw=itw+1; 
		
    if(itw>(log(optim.TolFun)/log(optim.tau3))),
      % No new optium found, linesearch failed.
      data.bracket_exitflag=-2; break; 
    end
    
    if(data.beta>0&&hill)
            % Get the brackets around minimum point
            % Pick bracket A from stored trials
            [t,i]=sort(data.storex,'ascend');
            storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
            [t,i]=find(storex>data.beta,1);
            if(isempty(i)), [t,i]=find(storex==data.beta,1); end
            alpha=storex(i); f_alpha=storefx(i); fPrime_alpha=storepx(i);
            
            % Pick bracket B from stored trials
            [t,i]=sort(data.storex,'descend');
            storefx=data.storefx(i);storepx=data.storepx(i); storex=data.storex(i);
            [t,i]=find(storex<data.beta,1);
            if(isempty(i)), [t,i]=find(storex==data.beta,1); end
            beta=storex(i); f_beta=storefx(i); fPrime_beta=storepx(i);
            
            % Calculate derivatives if not already calculated
            if(optim.GradConstr)
                gstep=data.initialStepLength/1e6; 
                if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
                if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
                [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
                [data,f_beta2]=gradient_function(data.xInitial(:)+(beta+gstep)*data.dir(:),funfcn, data, optim);
                fPrime_alpha=(f_alpha2-f_alpha)/gstep;
                fPrime_beta=(f_beta2-f_beta)/gstep;
            end
            % Set the brackets A and B
            data.a=alpha; data.f_a=f_alpha; data.fPrime_a=fPrime_alpha;
            data.b=beta; data.f_b=f_beta; data.fPrime_b=fPrime_beta;
  
            % Finished bracketing phase
            data.bracket_exitflag  = 2; return
    end
	% Reached max function evaluations
	if(data.funcCount>=optim.MaxFunEvals), data.bracket_exitflag=0; return; end
end
