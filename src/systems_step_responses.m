function systems_step_responses(systems)

    fns = fieldnames(systems);

    sys0        = systems.(fns{1});
    sysh_ident  = systems.(fns{2});
    sysh_kung   =  systems.(fns{3});
    
    figure
    t = 1:1:100;
    step(sys0,t)
    hold on
    step(sysh_ident,t)
    step(sysh_kung,t)
    legend('real system', 'Ident system', 'Kung system')

end