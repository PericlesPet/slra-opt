function M = sys_comparison(u, y, sys_ID, plot, t_id)

if exist('plot')
    if plot
        figure
        compare(iddata(y, u), idss(sys_ID)); 
    end
end

[Yh, M] = compare(iddata(y, u), idss(sys_ID)); 


if nargin>4
    X = sprintf('Time: %f\n', t_id);
    disp(X)
end

fprintf('\nAvg. Acc:  %3.3f%%\n', mean(M));
fprintf('Accuracy: [%3.3f%%]\n', M);


end