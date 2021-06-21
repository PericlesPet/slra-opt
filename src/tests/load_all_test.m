
% Load ALL Data
all_files = dir(fullfile('data\final Data\','*.mat'))
for i = 1:numel(all_files)
    file_path = [all_files(i).folder '\' all_files(i).name];
    fprintf('Loading: %s \n',file_path);
    load(file_path)
end