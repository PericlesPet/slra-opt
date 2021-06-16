function message=getexitmessage(exitflag)
    switch(exitflag)
        case 1, message='Change in the objective function value was less than the specified tolerance TolFun.';
        case 2, message='Change in x was smaller than the specified tolerance TolX.'; 
        case 3, message='Magnitude of gradient smaller than the specified tolerance';
        case 4, message='Boundary fminimum reached.';
        case 0, message='Number of iterations exceeded options.MaxIter or number of function evaluations exceeded options.FunEvals.';
        case -1, message='Algorithm was terminated by the output function.';
        case -2, message='Line search cannot find an acceptable point along the current search';
        otherwise, message='Undefined exit code';
    end
 