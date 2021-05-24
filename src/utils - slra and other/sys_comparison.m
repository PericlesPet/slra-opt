function M = sys_comparison(u, y, sys_ID, t_id, noPlot)

if exist('noPlot')
    if noPlot
        [Yh, M] = compare(iddata(y, u), idss(sys_ID)); 
    end
else
    figure
    compare(iddata(y, u), idss(sys_ID)); 
    [Yh, M] = compare(iddata(y, u), idss(sys_ID)); 
end

if nargin>3
    X = sprintf('Time: %f\n', t_id);
    disp(X)
end

X = sprintf('Accuracy : %f\n', M);
disp(X)


end