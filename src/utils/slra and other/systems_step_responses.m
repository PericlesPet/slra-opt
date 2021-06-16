function systems_step_responses(systems)

    fns = fieldnames(systems);

    sys0        = systems.sys0;
    sysh_ident  = systems.sysh_ident;
    if exist( 'systems.sysh_kung'), sysh_kung   =  systems.sysh_kung;end
    
    figure
    t = 1:1:100;
    step(sys0,t)
    hold on
    step(sysh_ident,t)
%     step(sysh_kung,t)
    legend('real system', 'Ident system')

end