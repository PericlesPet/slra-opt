function data = sectioningPhase_simple(funfcn, data, optim)
% Get the brackets
brcktEndpntA=data.a; brcktEndpntB=data.b;
% Calculate minimum between brackets
[alpha,f_alpha_estimated] = pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,data.a,data.b,data.f_a,data.fPrime_a,data.f_b,data.fPrime_b,optim);  
if(isfield(data,'beta')&&(data.f_beta<f_alpha_estimated)), alpha=data.beta; end
[t,i]=find(data.storex==alpha,1);
if((~isempty(i))&&(~isnan(data.storegx(i))))
    f_alpha=data.storefx(i); grad=data.storegx(:,i);
else
    % Calculate the error and gradient for the next minimizer itteration
    [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
    if(isfield(data,'beta')&&(data.f_beta<f_alpha)), 
        alpha=data.beta; 
        if((~isempty(i))&&(~isnan(data.storegx(i))))
            f_alpha=data.storefx(i); grad=data.storegx(:,i);
        else
            [data,f_alpha, grad]=gradient_function(data.xInitial(:)+alpha*data.dir(:),funfcn, data,optim);
        end
    end
end
% Store values linesearch
data.storefx=[data.storefx f_alpha]; data.storex=[data.storex alpha];
fPrime_alpha = grad'*data.dir(:);
data.alpha=alpha; 
data.fPrime_alpha= fPrime_alpha; 
data.f_alpha= f_alpha;
data.grad=grad;
% Set the exit flag to succes   
data.section_exitflag=[];
