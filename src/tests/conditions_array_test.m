array = [1 1 2 1 1 3];
target_array = [1 2 3 3 4 5 6 6 6];

total_elems = sum(array);
test_plot_array = 1:total_elems;

my_array = zeros(1, total_elems);
j = 1;
tmp_array = array;

for i=1:length(array)
    while tmp_array(i) > 0
        my_array(j) = i;
        j = j + 1;
        tmp_array(i) = tmp_array(i) - 1;
    end
end

% my_array - target_array
%%
plot(my_array, test_plot_array)
xlim([-1 10])
ylim([0 10])