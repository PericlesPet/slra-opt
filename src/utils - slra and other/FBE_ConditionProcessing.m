function [my_array, fbe_low_array, fbe_high_array, fbe_higher_array] = ...
    FBE_ConditionProcessing(condition_array, fbe_low_index, fbe_high_index, fbe_higher_index)

%%
% VARS: 
% condition_array
% iter_steps, bcktracking_steps, 

bcktracking_steps   = size(condition_array , 1);
iter_steps          = size(condition_array , 2);
num_vars            = size(condition_array , 3);

num_nonzero_array = zeros(iter_steps,1);
a = condition_array;
for i = 1:iter_steps
    num_nonzero_array(i) = numberofelements(nonzeros(a(:,i,1)));
end
%%
% e.g.
% array = [1 1 2 1 1 3];
% target_array = [1 2 3 3 4 5 6 6 6];
array       = num_nonzero_array;
tmp_array   = array;
total_elems = sum(array);
my_array = zeros(1, total_elems);
j = 1;

for i=1:length(array)
    while tmp_array(i) > 0
        my_array(j) = i;
        j = j + 1;
        tmp_array(i) = tmp_array(i) - 1;
    end
end

%%
slices = zeros(total_elems, 6);
slice = zeros(iter_steps*bcktracking_steps, 6);
for i = 1:4 %
    slice(:,i) = reshape(condition_array(:,:, i), iter_steps*bcktracking_steps, 1);
    slices(:,i) = nonzeros(slice(:, i));
end

%%
% fbe_low_index         = 1;
% fbe_high_index        = 4;
fbe_low_array   = slices(:,fbe_low_index) ;
fbe_high_array  = slices(:,fbe_high_index) ;

if nargin > 3 && nargout > 3
    fbe_higher_array  = slices(:,fbe_higher_index) ;
end

% plot(my_array, slices(:,fbe_low_index))
% hold on 
% plot(my_array, slices(:,fbe_high_index))

% 
% for subplot_4 = 1
%     plot(f_evals, 'ko', 'MarkerSize', 2)
%     hold on
%     plot(condition_array(1,:),'b')
%     plot(condition_array(4,:),'r')
%     legend('f evals', 'fbe 1', 'fbe 2')
% end
end
