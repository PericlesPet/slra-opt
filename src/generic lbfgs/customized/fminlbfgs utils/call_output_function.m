   
function stopt=call_output_function(data,optim,where)
stopt=false;
if(~isempty(optim.OutputFcn))
    output.iteration = data.iteration;
    output.funccount = data.funcCount;
    output.fval = data.fInitial;
    output.stepsize = data.alpha;
    output.directionalderivative = data.fPrimeInitial;
    output.gradient = reshape(data.gradient, data.xsizes);
    output.searchdirection = data.dir;
    stopt=feval(optim.OutputFcn,reshape(data.xInitial,data.xsizes),output,where); 
end
