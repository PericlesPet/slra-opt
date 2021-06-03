
% selectVisualize:
% 1 :   Gradient Descent Data
% 2 :   myFmincon        Data
% 3 :   fminlbfgs        Data 
selectVisualize = 3;

switch selectVisualize

    case 1
        
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


    case 2
        
    %% FROM myFMINCON
    % %% VISUALIZE FVALS, SLRA, CONSTRAINTS, AND ACCURACY. ALSO COMPARE TO OPTIMAL SLRA_MEX VALUES
    % ----> dataVisualize 2ND PART

    % fminconData = fminconData_aLM;
    % fminconData = fminconData_gdTrue;
    % fminconData = fminconData_gdfalse;
    fminconData = fminconData_slraVSslra
    use_iter = 0;
    if ~isfield(fminconData, 't_stamps'), ...
            fminconData.t_stamps = 1:length(fminconData.f_slra_val); use_iter = 1;end


    set(0,'DefaultLegendAutoUpdate','off')
    figure
    subplot(2,2,1)
    plot(fminconData.t_stamps, fminconData.f_slra_val)
    title('SLRA M(R)')
    hold on
    if ~use_iter, plot(info.iterinfo(1,:), info.iterinfo(2,:), 'r--'), else ...
        plot(1:length(info.iterinfo(1,:)), info.iterinfo(2,:), 'r--'), end
    line([0 fminconData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')
    legend('fmincon', 'slra mex', 'optimum point')

    subplot(2,2,2)
    for ii = 1:size(fminconData.CE, 1)    
        semilogy(fminconData.t_stamps, fminconData.CE(ii,:))
        hold on
    end
    % hold on
    % semilogy(fminconData.CE(2,:))
    title('Constraints')
    if size(fminconData.CE, 1) == 1
        legend('RR^T - I_N = 0 Const.')
    else
        legend('R*Hank(p_{hat}) = 0 Const.', 'RR^T - I_N = 0 Const.')
    end

    subplot(2,2,3)
    plot(fminconData.t_stamps, fminconData.fvals)
    title('Fvals')

    subplot(2,2,4)
    plot(fminconData.t_stamps, max(fminconData.M0, -100), 'b')
    % ylim([-200 100])
    title('Accuracy (%)')
    if size(fminconData.M0, 1) == 2
        line([0 fminconData.t_stamps(end)], ...
            [M_slra(1) M_slra(1)], 'Color', 'r')
        line([0 fminconData.t_stamps(end)], ...
            [M_slra(2) M_slra(2)], 'Color', 'r')
        legend('accuracy for y1', 'accuracy for y2', 'slra best accuracy for y1','slra best accuracy for y2')
    else
        line([0 fminconData.t_stamps(end)], ...
            [mean(M_slra) mean(M_slra)], 'Color', 'r')
        legend('accuracy for fmincon', 'slra best accuracy')
    end

    % subplot(2,2,3)
    % plot(fminconData.f_slra_val)

    %% VISUALIZE DIFFERENCE BETWEEN FVALS, SLRA, AND CONSTRAINTS (PROBABLY USELESS)
    for visualizeDiff = 1
        figure
        %%%%% SLRA %%%%%
        subplot(2,2,1)
        plot(fminconData.t_stamps, fminconData.f_slra_val)
        title('SLRA M(R)')
        hold on
        if ~use_iter, plot(info.iterinfo(1,:), info.iterinfo(2,:), 'r--'), else ...
            plot(1:length(info.iterinfo(1,:)), info.iterinfo(2,:), 'r--'), end
        line([0 fminconData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')
        legend('fmincon', 'slra mex', 'optimum point')

        %%%%% FVALS %%%%%
        subplot(2,2,3)
        plot(fminconData.t_stamps, 2*fminconData.fvals)
        title('Fvals')

        %%%%% SLRA - FVALS %%%%%
        subplot(2,2,2)
        plot(fminconData.t_stamps, fminconData.f_slra_val-2*fminconData.fvals)
        title('SLRA - Fvals')

        %%%%% CONSTRAINTS%%%%%
        subplot(2,2,4)
        for ii = 1:size(fminconData.CE, 1)    
            semilogy(fminconData.t_stamps, fminconData.CE(ii,:))
            hold on
        end
        % hold on
        % semilogy(fminconData.CE(2,:))
        title('Constraints')
        if size(fminconData.CE, 1) == 1
            legend('RR^T - I_N = 0 Const.')
        else
            legend('R*Hank(p_{hat}) = 0 Const.', 'RR^T - I_N = 0 Const.')
        end
    end

    
    
    case 3
    %% VISUALIZE FMINLBFGS 
    % (difference: 2 Xs / f_vals for lbfgs and prox 
    % --> data.fvals_fminlbfgs      / Xs_lbfgs
    % --> data.fvals_fminlbfgsProx  / Xs_prox
    
    % %% VISUALIZE FVALS, SLRA, CONSTRAINTS, AND ACCURACY. ALSO COMPARE TO OPTIMAL SLRA_MEX VALUES 
    % ----> dataVisualize 2ND PART

        % fminconData = fminconData_aLM;
        % fminconData = fminconData_gdTrue;
%         fminconData = fminconData_gdfalse;
        lbfgsData = fminlbfgsData
        use_iter = 0;
        if ~isfield(lbfgsData, 't_stamps'), ...
                lbfgsData.t_stamps = 1:length(lbfgsData.f_slra_val); use_iter = 1;end


        set(0,'DefaultLegendAutoUpdate','off')
        figure
        subplot(2,2,1)
        plot(lbfgsData.t_stamps, lbfgsData.f_slra_val)
        title('SLRA M(R)')
        hold on
        if ~use_iter, plot(info.iterinfo(1,:), info.iterinfo(2,:), 'r--'), else ...
            plot(1:length(info.iterinfo(1,:)), info.iterinfo(2,:), 'r--'), end
        line([0 lbfgsData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')
        legend('fmincon', 'slra mex', 'optimal slra')
        ylabel('M(R)')
        xlabel('t (seconds)')

        subplot(2,2,2)
        for ii = 1:size(lbfgsData.CE, 1)    
            semilogy(lbfgsData.t_stamps, lbfgsData.CE(ii,:))
            hold on
        end
        % hold on
        % semilogy(lbfgsData.CE(2,:))
        title('Constraints')
        if size(lbfgsData.CE, 1) == 1
            legend('RR^T - I_N = 0 Const.')
        else
            legend('R*Hank(p_{hat}) = 0 Const.', 'RR^T - I_N = 0 Const.')
        end
        xlabel('t (seconds)')

        subplot(2,2,3)
        plot(lbfgsData.t_stamps, lbfgsData.fvals)
        title('Fvals')
        ylabel('f(x)')
        xlabel('t (seconds)')

        subplot(2,2,4)
        plot(lbfgsData.t_stamps, max(lbfgsData.M0, -100), 'b')
        % ylim([-200 100])
        title('Accuracy (%)')
        if size(lbfgsData.M0, 1) == 2
            line([0 lbfgsData.t_stamps(end)], ...
                [M_slra(1) M_slra(1)], 'Color', 'r')
            line([0 lbfgsData.t_stamps(end)], ...
                [M_slra(2) M_slra(2)], 'Color', 'r')
            legend('accuracy for y1', 'accuracy for y2', 'slra best accuracy for y1','slra best accuracy for y2')
        else
            line([0 lbfgsData.t_stamps(end)], ...
                [mean(M_slra) mean(M_slra)], 'Color', 'r')
            legend('accuracy for fmincon', 'slra best accuracy')
        end
        ylabel('accuracy (%)')
        xlabel('t (seconds)')

        % subplot(2,2,3)
        % plot(lbfgsData.f_slra_val)
end


        
        
        