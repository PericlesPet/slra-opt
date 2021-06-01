

% Proj Data
max(mean(gdProjData.M(:,:)));
semilogy(gdProjData.residual_R)
plot(gdProjData.t_stamps, gdProjData.f)
plot(gdProjData.t_stamps, gdProjData.M)
line([t_slra t_slra], [86 100], 'Color', 'r')
line([0 7], [gdProjData.Mopt(1) gdProjData.Mopt(1)], 'Color', 'r')
line([0 7], [gdProjData.Mopt(2) gdProjData.Mopt(2)], 'Color', 'r')
plot(info.iterinfo(1,:), info.iterinfo(2,:))
hold on 
plot(gdProjData.t_stamps, gdProjData.f)



% Reg Data
plot(info.iterinfo(1,:), info.iterinfo(2,:))
hold on 
plot(gdRegData.t_stamps, gdRegData.f)
