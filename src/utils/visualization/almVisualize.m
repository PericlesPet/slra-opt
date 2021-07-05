%% Extract data into vectors
plotALMdata = almData;
fvals           = plotALMdata.fval;
xs              = plotALMdata.x;
lambdas         = plotALMdata.lambda;
vs              = plotALMdata.v;
rhos            = plotALMdata.rho;
kkt_1           = [plotALMdata.kkt{1,:}];
kkt_2           = [plotALMdata.kkt{2,:}];
kkt_3           = [plotALMdata.kkt{3,:}];

Ls              = [plotALMdata.searchData(:).L];
Mslras          = [plotALMdata.searchData(:).Mslra];
CEs             = [plotALMdata.searchData(:).CE];
DPs             = [plotALMdata.searchData(:).DP];
% DPus            = [checkdata.searchData(:).DPu];
Ms              = [plotALMdata.searchData(:).M];
iters           = [plotALMdata.searchData(:).iters];
t_stamps        = [plotALMdata.searchData(:).t_stamps];

figure 
subplot(3,2,1)
plot(t_stamps,Ls)
xlabel('t (sec)')
title('L(x,R,y)')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.9* min(Ls) 1.1*max(Ls)], 'LineStyle' , ':', 'Color', 'r')
end    

subplot(3,2,3)
semilogy(t_stamps,CEs')
xlabel('t (sec)')
title('norm(RH(x))     &     norm(R*R^{T} - I)')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.9* min(min(CEs)) 1.1*max(max(CEs))], 'LineStyle' , ':', 'Color', 'r')
end    
legend('norm(RH(x))', 'norm(R*R^{T} - I)')

subplot(3,2,5)
plot(t_stamps,DPs)
xlabel('t (sec)')
title('norm(p - p_{hat})')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.9* min(DPs) 1.1*max(DPs)], 'LineStyle' , ':', 'Color', 'r')
end

subplot(3,2,2)
plot(t_stamps,Mslras)
xlabel('t (sec)')
title('M(R) (Blue) Vs. Optimal M(R*) (Red)')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.8* min(Mslras) 1.1*max(Mslras)], 'LineStyle' , ':', 'Color', 'r')
end    
line([0 t_stamps(x_curr)], [info_mex.fmin info_mex.fmin], 'Color', 'r')
ylim([0.8* min(Mslras) 1.1*max(Mslras)])

subplot(3,2,4)

plot(t_stamps,mean(Ms))
xlabel('t (sec)')
title('Average Accuracy (Blue) Vs. Optimal Accuracy (Red) (%) ')
x_curr = 0;
for ii = 1:length(iters)
	x_curr = x_curr + iters(ii);
    line([t_stamps(x_curr) t_stamps(x_curr)], [0.9* min(mean(Ms)) 1.1*max(mean(Ms))], 'LineStyle' , ':', 'Color', 'r')
end
[~, M1_opt] = sysAccuracy(info_mex.Rh);
line([0 t_stamps(x_curr)], [mean(M1_opt) mean(M1_opt)], 'Color', 'r')
ylim([80 100])


suptitle({'ALM','(vertical lines mark outer iterations)'})

% subplot(3,2,6)
% plot(DPus)
% title('norm(u - u_{hat})')
% x_curr = 0;
% for ii = 1:length(iters)
% 	x_curr = x_curr + iters(ii);
%     line([x_curr x_curr], [0.9* min(DPus) 1.1*max(DPus)], 'LineStyle' , ':', 'Color', 'r')
% end


% legend('1', '2', '3', '4','5')

