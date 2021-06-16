function data = bracketingPhase(funfcn, data, optim)
% bracketingPhase finds a bracket [a,b] that contains acceptable points; a bracket 
% is the same as a closed interval, except that a > b is allowed.
%
% The outputs f_a and fPrime_a are the values of the function and the derivative 
% evaluated at the bracket endpoint 'a'. Similar notation applies to the endpoint 
% 'b'. 
% Parameters of bracket A
data.a = []; 
data.f_a = []; 
data.fPrime_a = []; 
% Parameters of bracket B
data.b = []; 
data.f_b = []; 
data.fPrime_b = [];
% First trial alpha is user-supplied
% f_alpha will contain f(alpha) for all trial points alpha
% fPrime_alpha will contain f'(alpha) for all trial points alpha
alpha = data.initialStepLength;
f_alpha = data.fInitial;              
fPrime_alpha = data.fPrimeInitial;    
% Set maximum value of alpha (determined by fminimum)
alphaMax = (data.fminimum - data.fInitial)/(optim.rho*data.fPrimeInitial); 
alphaPrev = 0;
while(true) 
  % Evaluate f(alpha) and f'(alpha)
  fPrev = f_alpha;
  fPrimePrev = fPrime_alpha;
  
  % Calculate value (and gradient if no extra time cost) of current alpha
  if(~optim.GradConstr)
      [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
      fPrime_alpha = grad'*data.dir(:);
  else
      gstep=data.initialStepLength/1e6;
      if(gstep>optim.DiffMaxChange), gstep=optim.DiffMaxChange; end
      if(gstep<optim.DiffMinChange), gstep=optim.DiffMinChange; end
      [data,f_alpha]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
      [data,f_alpha2]=gradient_function(data.xInitial(:)+(alpha+gstep)*data.dir(:),funfcn, data, optim);
      fPrime_alpha=(f_alpha2-f_alpha)/gstep;
  end
  
  % Store values linesearch 
  data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha]; 
	
  % Terminate if f < fminimum
  if (f_alpha <= data.fminimum), data.bracket_exitflag = 4; return; end
  
  % Bracket located - case 1 (Wolfe conditions)
  if (f_alpha > (data.fInitial + alpha*optim.rho*data.fPrimeInitial)) || (f_alpha >= fPrev)
    % Set the bracket values
    data.a = alphaPrev; data.f_a = fPrev;  data.fPrime_a = fPrimePrev;
    data.b = alpha; data.f_b = f_alpha;  data.fPrime_b = fPrime_alpha;
    % Finished bracketing phase
    data.bracket_exitflag  = 2; return 
  end
  % Acceptable steplength found
  if (abs(fPrime_alpha) <= -optim.sigma*data.fPrimeInitial), 
      if(optim.GradConstr)
          % Gradient was not yet calculated because of time costs
          [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data, optim);
          fPrime_alpha = grad'*data.dir(:);
      end
      % Store the found alpha values
      data.alpha=alpha;
      data.fPrime_alpha= fPrime_alpha; data.f_alpha= f_alpha; data.grad=grad;
      % Finished bracketing phase, and no need to call sectioning phase
      data.bracket_exitflag = [];  return 
  end
  
  % Bracket located - case 2  
  if (fPrime_alpha >= 0)
    % Set the bracket values
    data.a = alpha; data.f_a = f_alpha;  data.fPrime_a = fPrime_alpha;
    data.b = alphaPrev; data.f_b = fPrev; data.fPrime_b = fPrimePrev;
    % Finished bracketing phase
    data.bracket_exitflag  = 2; return
  end
 
  % Update alpha
  if (2*alpha - alphaPrev < alphaMax )
      brcktEndpntA = 2*alpha-alphaPrev; 
      brcktEndpntB = min(alphaMax,alpha+optim.tau1*(alpha-alphaPrev));
      % Find global minimizer in bracket [brcktEndpntA,brcktEndpntB] of 3rd-degree polynomial 
      % that interpolates f() and f'() at alphaPrev and at alpha
      alphaNew = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alphaPrev,alpha,fPrev, ...
                                         fPrimePrev,f_alpha,fPrime_alpha,optim);
      alphaPrev = alpha;
      alpha = alphaNew;
  else
      alpha = alphaMax;
  end
  % maxFunEvals reached
  if(data.funcCount >optim.MaxFunEvals), data.bracket_exitflag = -1; return, end
end
