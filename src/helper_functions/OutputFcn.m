function stop = OutputFcn(x, optimValues, state) 
global Th F; Th = [Th x']; F = [F optimValues.fval']; stop = 0;
