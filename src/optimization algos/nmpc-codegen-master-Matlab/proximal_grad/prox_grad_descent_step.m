function [ x_new,gamma ] = prox_grad_descent_step( x,gamma,beta,proxg,f,df )    
    x_new = proxg(x-gamma*df(x));
    r=x_new-x;
%     
%     value1 = f(x_new) 
% %     value2_all = f(x) - dot(df(x),x_new-x) + norm(x_new-x,2)
%     value2_all = f(x) + dot(df(x),x_new-x) + norm(x_new-x,2)
%     value2 = f(x)
%     value3 = - dot(df(x),x_new-x)
% %     value3 = +(1-beta)/(2*gamma) * norm(x_new-x,2)
%     value4 = norm(x_new-x,2)
    while(f(x_new)>f(x)- dot(df(x),x_new-x)+(1-beta)/(2*gamma) * norm(x_new-x,2))
        gamma=gamma/2;
        x_new = proxg(x-gamma*df(x));
    end
%     value1 = f(x_new) 
%     value2 = f(x) - dot(df(x),x_new-x) + norm(x_new-x,2)
%     while(f(x_new)>f(x) + dot(df(x),x_new-x) + norm(x_new-x,2))
%     while(f(x_new)>f(x))
%         inner_value1 = f(x_new) 
%         inner_value2 = f(x) + dot(df(x),x_new-x) + norm(x_new-x,2)
%         gamma=gamma/2
%         x_new = proxg(x-gamma*df(x));
%     end
end