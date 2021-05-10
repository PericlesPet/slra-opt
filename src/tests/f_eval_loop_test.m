
f_evals = zeros(number_of_steps,1);
for i = 1:number_of_steps
    f_evals(i) = f(x_steps(:,i));
end
plot( 1:number_of_steps, f_evals)
    
    