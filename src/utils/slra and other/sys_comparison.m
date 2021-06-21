function M = sys_comparison(u, y, sys_ID, plot, t_id)

if exist('plot')
    if plot
        figure
        compare(iddata(y, u), idss(sys_ID)); 
    end
end

[Yh, M] = compare(iddata(y, u), idss(sys_ID)); 
M = max(M, -100);


fprintf('Avg. Acc:  %3.3f%%\n', mean(M));
fprintf('Accuracy: [%3.3f%%]\n', M);
if nargin>4
    fprintf('Sys Accuracy Time: %f\n', t_id);
end
fprintf('\n');


end