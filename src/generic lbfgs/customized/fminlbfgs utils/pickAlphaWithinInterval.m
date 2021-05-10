function [alpha,f_alpha]= pickAlphaWithinInterval(brcktEndpntA,brcktEndpntB,alpha1,alpha2,f1,fPrime1,f2,fPrime2,optim)
% finds a global minimizer alpha within the bracket [brcktEndpntA,brcktEndpntB] of the cubic polynomial 
% that interpolates f() and f'() at alpha1 and alpha2. Here f(alpha1) = f1, f'(alpha1) = fPrime1, 
% f(alpha2) = f2, f'(alpha2) = fPrime2.
% determines the coefficients of the cubic polynomial with c(alpha1) = f1, 
% c'(alpha1) = fPrime1, c(alpha2) = f2, c'(alpha2) = fPrime2.
coeff = [(fPrime1+fPrime2)*(alpha2-alpha1)-2*(f2-f1) ...
    3*(f2-f1)-(2*fPrime1+fPrime2)*(alpha2-alpha1) (alpha2-alpha1)*fPrime1 f1];
% Convert bounds to the z-space
lowerBound = (brcktEndpntA - alpha1)/(alpha2 - alpha1);
upperBound = (brcktEndpntB - alpha1)/(alpha2 - alpha1);
% Swap if lowerbound is higher than the upperbound
if (lowerBound  > upperBound), t=upperBound; upperBound=lowerBound; lowerBound=t; end 
% Find minima and maxima from the roots of the derivative of the polynomial.
sPoints = roots([3*coeff(1) 2*coeff(2) coeff(3)]); 
% Remove imaginaire and points outside range
sPoints(imag(sPoints)~=0)=[]; 
sPoints(sPoints<lowerBound)=[]; sPoints(sPoints>upperBound)=[];
% Make vector with all possible solutions
sPoints=[lowerBound sPoints(:)' upperBound];
% Select the global minimum point
[f_alpha,index]=min(polyval(coeff,sPoints)); z=sPoints(index);
% Add the offset and scale back from [0..1] to the alpha domain
alpha = alpha1 + z*(alpha2 - alpha1);
% Show polynomial search
if(optim.Display(1)=='p'); 
    vPoints=polyval(coeff,sPoints);
    plot(sPoints*(alpha2 - alpha1)+alpha1,vPoints,'co');
    plot([sPoints(1) sPoints(end)]*(alpha2 - alpha1)+alpha1,[vPoints(1) vPoints(end)],'c*');
    xPoints=linspace(lowerBound/3, upperBound*1.3, 50);
    vPoints=polyval(coeff,xPoints);
    plot(xPoints*(alpha2 - alpha1)+alpha1,vPoints,'c');
end
