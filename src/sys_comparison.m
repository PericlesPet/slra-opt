function sys_comparison(u, y, sys_ID, t_id)

figure
compare(iddata(y, u), idss(sys_ID)); 
[Yh, M] = compare(iddata(y, u), idss(sys_ID)); 


X = sprintf('Time: %f\n', t_id);
disp(X)
X = sprintf('Accuracy : %f\n', M);
disp(X)


end