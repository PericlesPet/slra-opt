% complexitiesVisualize is for any visualizations to do with 
% Problem Complexity and Scaling

% Not to be confused with dataVisualizeComplexities, which is for data visualization while running 
% the algorithms for various complexties
 
%% Concatenate all tables/experiments
Alldata = [statsTable0, statsTable1, statsTable2];


concatenateddata = struct;  %final structure
for field = fieldnames(Alldata)'
   fname = field{1};
%    concatenateddata.(fname) = vertcat(Alldata.(fname));
   concatenateddata.(fname) = horzcat(Alldata.(fname));
end
statsTable = concatenateddata;


%%   sort data
slratimes = [statsTable.slra.t_ident];
[~, sortedI] = sort(slratimes)

sortedData = struct;  %final structure
for field = fieldnames(statsTable)'
   fname = field{1};
%    concatenateddata.(fname) = vertcat(Alldata.(fname));
   sortedData.(fname) = statsTable.(fname)(sortedI);
end


%% extract Ms, Ts, 

sortedData = struct;  %final structure
for field = fieldnames(statsTable)'
   fname = field{1};
%    concatenateddata.(fname) = vertcat(Alldata.(fname));
   sortedData.(fname) = statsTable.(fname)(sortedI);
end

