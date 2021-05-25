%% Extract data into vectors
fvals           = checkdata.fval;
xs              = checkdata.x;
lambdas         = checkdata.lambda;
vs              = checkdata.v;
rhos            = checkdata.rho;
kkt_1           = [checkdata.kkt{1,:}];
kkt_2           = [checkdata.kkt{2,:}];
kkt_3           = [checkdata.kkt{3,:}];

Ls              = [checkdata.searchData(:).L];
Mslras          = [checkdata.searchData(:).Mslra];
CEs             = [checkdata.searchData(:).CE];
DPs             = [checkdata.searchData(:).DP];
Ms              = [checkdata.searchData(:).M];
iters           = [checkdata.searchData(:).iters];


subplot(3,2,1)
plot(Ls)
title('L(x,R,y)')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([x_curr x_curr], [0.9* min(Ls) 1.1*max(Ls)], 'LineStyle' , ':', 'Color', 'r')
end    

subplot(3,2,3)
semilogy(CEs)
title('norm(RH(x))')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([x_curr x_curr], [0.9* min(CEs) 1.1*max(CEs)], 'LineStyle' , ':', 'Color', 'r')
end    


subplot(3,2,5)
plot(DPs)
title('norm(p - p_{hat})')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([x_curr x_curr], [0.9* min(DPs) 1.1*max(DPs)], 'LineStyle' , ':', 'Color', 'r')
end

subplot(3,2,2)
plot(Mslras)
title('M(R) (Blue) Vs. Optimal M(R*) (Red)')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([x_curr x_curr], [0.9* min(Mslras) 1.1*max(Mslras)], 'LineStyle' , ':', 'Color', 'r')
end    
line([0 x_curr], [info1.fmin info1.fmin], 'Color', 'r')


subplot(3,2,4)

plot(mean(Ms))
title('Average Accuracy (Blue) Vs. Optimal Accuracy (Red) (%) ')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([x_curr x_curr], [0.9* min(mean(Ms)) 1.1*max(mean(Ms))], 'LineStyle' , ':', 'Color', 'r')
end
[~, M1_opt] = accuracy_r(info1.Rh);
line([0 x_curr], [mean(M1_opt) mean(M1_opt)], 'Color', 'r')

% legend('1', '2', '3', '4','5')

