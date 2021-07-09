% complexitiesVisualize is for any visualizations to do with 
% Problem Complexity and Scaling

% Not to be confused with dataVisualizeComplexities, which is for data visualization while running 
% the algorithms for various complexties

%% Load StatsTable
all_files = dir(fullfile('data\final Data\','*.mat'));
for i = 1:numel(all_files)
    if strcmp(all_files(i).name(1:end-4), 'statsTable')
        file_path = [all_files(i).folder '\' all_files(i).name];
        fprintf('Loading: %s \n',file_path);
        load(file_path)
    end
end



%%   sort data
slratimes = [statsTable.slra.t_ident];
[~, sortedI] = sort(slratimes);

sortedData = struct;  %final structure
for field = fieldnames(statsTable)'
   fname = field{1};
%    concatenateddata.(fname) = vertcat(Alldata.(fname));
   sortedData.(fname) = statsTable.(fname)(sortedI);
end


% extract Ms, Ts, 

sortedData = struct;  %final structure
for field = fieldnames(statsTable)'
   fname = field{1};
%    concatenateddata.(fname) = vertcat(Alldata.(fname));
   sortedData.(fname) = statsTable.(fname)(sortedI);
end


%%
clc

for index = 1:10

    fprintf('\nSTATS FOR EXP %d\nCOMPLEXITY = %d \n\n', index, sortedData.complexities(index))
    slrat = sortedData.slra(index).t_ident
    % slraM = mean(sortedData.slra(index).    
    fminlbfgs_prox_exp   = sortedData.fminlbfgs_prox(index)
    fminlbfgs_simple_exp = sortedData.fminlbfgs_simple(index)
    panocfminlbfgs_exp   = sortedData.panocfminlbfgs(index)

    fmincon_exp     = sortedData.fmincon(index)
    gd_reg_exp      = sortedData.gd_reg(index)
    gd_proj_exp     = sortedData.gd_proj(index)
    panoclbfgs_exp  = sortedData.panoclbfgs(index)

end




%% Plot Times for various Complexities
figure
plot(ctimes.Fmincon, 'k-o')
hold on 
plot(ctimes.Gdreg, '--s')
plot(ctimes.Gdproj, '--d')
% plot(ctimes.Panoc, '-..')
plot(ctimes.Slra, 'r-+')
set(gca,'YScale','log')
legend('fmincon', ...
    'GD_{reg}', 'GD_{proj}', ...
    'SLRA C++')

ylabel('Logarithmic Time (seconds)') 
xlabel('Experiment Number')
title('SLRA C++ Vs Gradient Methods & FMINCON')


figure
plot(ctimes.Fminlbfgsprox, '--p')
hold on 
plot(ctimes.Fminlbfgssimple, '--^')
plot(ctimes.Panocfminlbfgs, '-.*')
plot(ctimes.Slra, 'r-+')
plot(ctimes.Panoc, '-..')

ylabel('Logarithmic Time (seconds)') 
xlabel('Experiment Number')
title('SLRA C++ Vs Proximal Methods')

set(gca,'YScale','log')
legend('fminlbfgs_{Prox}', 'fminlbfgs_{Plain}', ...
    'PANOC_{fminlbfgs}', 'SLRA C++', 'PANOC')
% sortedData.



% %% Concatenate all tables/experiments
% Alldata = [statsTable0, statsTable1, statsTable2];
% 
% 
% concatenateddata = struct;  %final structure
% for field = fieldnames(Alldata)'
%    fname = field{1};
% %    concatenateddata.(fname) = vertcat(Alldata.(fname));
%    concatenateddata.(fname) = horzcat(Alldata.(fname));
% end
% statsTable = concatenateddata;
