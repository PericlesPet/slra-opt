function data = sectioningPhase(funfcn, data, optim)
%
% sectioningPhase finds an acceptable point alpha within a given bracket [a,b] 
% containing acceptable points. Notice that funcCount counts the total number of 
% function evaluations including those of the bracketing phase. 
while(true)
    
    % Pick alpha in reduced bracket
    brcktEndpntA = data.a + min(optim.tau2,optim.sigma)*(data.b - data.a); 
    brcktEndpntB = data.b - optim.tau3*(data.b - data.a);
    
    % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree 
    % polynomial that interpolates f() and f'() at "a" and at "b".
    alpha = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);  
    % No acceptable point could be found
    if (abs( (alpha - data.a)*data.fPrime_a ) <= data.TolFunLnS), data.section_exitflag = -2; return; end
    
    % Calculate value (and gradient if no extra time cost) of current alpha
    if(~optim.GradConstr)
        [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
        fPrime_alpha = grad'*data.dir(:);
    else
        gstep=data.initialStepLength/1e6; 
        if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
        if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
        [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
        [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
        fPrime_alpha=(f_alpha2-f_alpha)/gstep;
    end
	% Store values linesearch 
	data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha]; 
	
    % Store current bracket position of A
    aPrev = data.a; 
    f_aPrev = data.f_a; 
    fPrime_aPrev = data.fPrime_a; 
    % Update the current brackets
    if ((f_alpha > data.fInitial + alpha*optim.rho*data.fPrimeInitial) || (f_alpha >= data.f_a))
        % Update bracket B to current alpha
        data.b = alpha; data.f_b = f_alpha; data.fPrime_b = fPrime_alpha;
    else
        % Wolfe conditions, if true then acceptable point found 
        if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial), 
            if(optim.GradConstr)
                % Gradient was not yet calculated because of time costs
                [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
                fPrime_alpha = grad'*data.dir(:);
            end
            % Store the found alpha values
            data.alpha=alpha; data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha;
            data.grad=grad;
            data.section_exitflag = []; return, 
        end
        
        % Update bracket A
        data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;
        
        if (data.b - data.a)*fPrime_alpha >= 0
            % B becomes old bracket A;
            data.b = aPrev; data.f_b = f_aPrev;  data.fPrime_b = fPrime_aPrev;
        end
    end
    
    % No acceptable point could be found
    if (abs(data.b-data.a) < eps), data.section_exitflag = -2; return, end
    % maxFunEvals reached
    if(data.funcCount >optim.MaxFunEvals), data.section_exitflag = -1; return, end
end
