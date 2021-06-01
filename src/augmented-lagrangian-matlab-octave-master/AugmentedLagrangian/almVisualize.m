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
% DPus            = [checkdata.searchData(:).DPu];
Ms              = [checkdata.searchData(:).M];
iters           = [checkdata.searchData(:).iters];
t_stamps        = [checkdata.searchData(:).t_stamps];

subplot(3,2,1)
plot(t_stamps,Ls)
title('L(x,R,y)')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.9* min(Ls) 1.1*max(Ls)], 'LineStyle' , ':', 'Color', 'r')
end    

subplot(3,2,3)
semilogy(t_stamps,CEs')
title('norm(RH(x))')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.9* min(min(CEs)) 1.1*max(max(CEs))], 'LineStyle' , ':', 'Color', 'r')
end    


subplot(3,2,5)
plot(t_stamps,DPs)
title('norm(p - p_{hat})')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.9* min(DPs) 1.1*max(DPs)], 'LineStyle' , ':', 'Color', 'r')
end

subplot(3,2,2)
plot(t_stamps,Mslras)
title('M(R) (Blue) Vs. Optimal M(R*) (Red)')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.9* min(Mslras) 1.1*max(Mslras)], 'LineStyle' , ':', 'Color', 'r')
end    
line([0 t_stamps(x_curr)], [info1.fmin info1.fmin], 'Color', 'r')
ylim([3 10])

subplot(3,2,4)

plot(t_stamps,mean(Ms))
title('Average Accuracy (Blue) Vs. Optimal Accuracy (Red) (%) ')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.9* min(mean(Ms)) 1.1*max(mean(Ms))], 'LineStyle' , ':', 'Color', 'r')
end
[~, M1_opt] = accuracy_r(info1.Rh);
line([0 t_stamps(x_curr)], [mean(M1_opt) mean(M1_opt)], 'Color', 'r')
ylim([80 100])


% subplot(3,2,6)
% plot(DPus)
% title('norm(u - u_{hat})')
% x_curr = 0;
% for ii = 1:length(iters)
% 	x_curr = x_curr + iters(ii);
%     line([x_curr x_curr], [0.9* min(DPus) 1.1*max(DPus)], 'LineStyle' , ':', 'Color', 'r')
% end


% legend('1', '2', '3', '4','5')

