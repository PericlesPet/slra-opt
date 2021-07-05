% Script for All Data Visualizations
% select Visual:
% 1 :   gradient descent Data
% 2 :   myFmincon        Data
% 3 :   fminlbfgs        Data 
% 4 :   panoc            Data
% 5 :   alm              Data

if ~exist('isAccSemilog'), isAccSemilog = 0; end

if isCloseAll == 1
    close all
    clc
end

% AGGREGATE FIGURES
aggrTstmps = [];
aggrMsMax  = [];
aggrMsMin  = [];

    % AGGREGATE ACCURACIES
aggrMfig = figure;
hold on

fgDetails = ['m_{in} = ', num2str(m_in), ', p_{out} = ',num2str(p_out), ...
    ', ell = ', num2str(ell) ,', T = ', num2str(T)];
title({'Aggregate Accuracies (%)', fgDetails})

    % AGGREGATE FVALS
aggrFVALfig = figure;
hold on
title({'Aggregate SLRA M(R) Vals', fgDetails})



for selectVisual = 1:4
    switch selectVisual

        case 1
        %% FROM GD
            gdXlim = 20;

            for selectGDesc = selectGDs
                if selectGDesc == 1
                    gdescData   = gdData;    
                    gdescTitle  = 'Simple Gradient Descent';
                    fgName      = 'gd_simple';
                    isAggrMfig    = 0;
                    isAggrFVALfig = 0;                    
                elseif selectGDesc == 4
                    gdescData   = gdRegData;    
                    gdescTitle  = 'Regularized Gradient Descent';
                    fgName      = 'gd_reg';
                    isAggrMfig    = 1;
                    isAggrFVALfig = 1;
                elseif selectGDesc == 2
                    gdescData   = gdManoptData;    
                    gdescTitle  = 'ManOpt-like Gradient Descent';
                    fgName      = 'gd_manopt';
                    isAggrMfig    = 0;
                    isAggrFVALfig = 0;
                elseif selectGDesc == 3
                    gdescData   = gdProjData;    
                    gdescTitle  = 'Stiefel-Projected Gradient Descent';
                    fgName      = 'gd_proj';
                    isAggrMfig    = 1;
                    isAggrFVALfig = 1;
                end

                use_iter = 0;
                if ~isfield(gdescData, 't_stamps'), ...
                        gdescData.t_stamps = 1:length(gdescData.f); use_iter = 1;end
                set(0,'DefaultLegendAutoUpdate','off')

                
                % AGGREGATE FIGURES 
                if isAggrMfig
                    aggrTstmps = [aggrTstmps gdescData.t_stamps(end)];
                    figure(aggrMfig)
                    if isAccSemilog
                        logMarray = max(mean(gdescData.M(:,:),1), ...
                            1./(abs(mean(gdescData.M(:,:),1))+1));
                        aggrMsMin = [aggrMsMin min(logMarray)];
                        aggrMsMax = [aggrMsMax max(logMarray)];                        
                        semilogy(gdescData.t_stamps, ...
                            logMarray)
                        ylim([min(aggrMsMin) 100])

                        %semilogy(gdescData.t_stamps, mean(gdescData.M))
                    else
                        plot(gdescData.t_stamps, max(mean(gdescData.M), 0))
                    end
                end
                if isAggrFVALfig
                    figure(aggrFVALfig)
                    plot(gdescData.t_stamps, gdescData.f)
                end
            end


        case 2
        %% FROM myFMINCON

        % %% VISUALIZE FVALS, SLRA, CONSTRAINTS, AND ACCURACY. ALSO COMPARE TO OPTIMAL SLRA_MEX VALUES
        % ----> dataVisualize 2ND PART

        for selectFminconFunc = selectFmincons
            if selectFminconFunc == 1
                fminconData = fminconData_aLM;
                fminconTitle = 'Matlab Fmincon for: L(x,R) s.t. RR^{T} - I = 0, compared to C SLRA';
                fgName      = 'fmincon_alm';
                isAggrMfig    = 0;
                isAggrFVALfig = 0;

            elseif selectFminconFunc == 3
                fminconData = fminconData_gdTrue;
                fminconTitle = 'Matlab Fmincon for: f(x) = |p-x|^{2} s.t. both constraints, compared to C SLRA';
                fgName      = 'fmincon_gdTrue';
                isAggrMfig    = 0;
                isAggrFVALfig = 0;
            elseif selectFminconFunc == 4
                fminconData = fminconData_gdFalse;
                fminconTitle = 'Matlab Fmincon for: f(x) = |p-x|^{2} s.t. both constraints, compared to C SLRA';
                fgName      = 'fmincon_gdFalse';
                isAggrMfig    = 0;
                isAggrFVALfig = 0;
            elseif selectFminconFunc == 2
                fminconData = fminconData_slraVSslra;
                fminconTitle = 'Matlab Fmincon for: M(R) s.t. RR^{T} - I = 0, compared to C SLRA';
                fgName      = 'fmincon_varpro';
                isAggrMfig    = 1;
                isAggrFVALfig = 1;
            end


            use_iter = 0;
            if ~isfield(fminconData, 't_stamps'), ...
                    fminconData.t_stamps = 1:length(fminconData.f_slra_val); use_iter = 1;end


            % AGGREGATE FIGURES 
            
            if isAggrMfig
                aggrTstmps = [aggrTstmps gdescData.t_stamps(end)];
                figure(aggrMfig)
                if isAccSemilog
                    logMarray = max(mean(fminconData.M0(:,:),1), ...
                            1./(abs(mean(fminconData.M0(:,:),1))+1));
                    aggrMsMin = [aggrMsMin min(logMarray)];
                    aggrMsMax = [aggrMsMax max(logMarray)];                        
                    semilogy(fminconData.t_stamps, ...
                        logMarray)
                    ylim([min(aggrMsMin) 100])
                    %semilogy(fminconData.t_stamps, mean(fminconData.M0))
                else                    
                    plot(fminconData.t_stamps, max(mean(fminconData.M0), 0))
                end
            end
            if isAggrFVALfig
                figure(aggrFVALfig)
                plot(fminconData.t_stamps, fminconData.f_slra_val)
            end
        end


        case 3
        %% FROM FMINLBFGS 

        % %% VISUALIZE FVALS, SLRA, CONSTRAINTS, AND ACCURACY. ALSO COMPARE TO OPTIMAL SLRA_MEX VALUES 
            for selectFminlbfgs = 1:2
                if selectFminlbfgs == 1
                    fminlbfgsData = fminlbfgsData_simple;
                    fminlbfgsTitle = 'fminlbfgs for M(R), compared to C SLRA';
                    fgName      = 'fminlbfgs_simple';

                elseif selectFminlbfgs == 2
                    fminlbfgsData = fminlbfgsData_prox;
                    fminlbfgsTitle = 'fminlbfgs for M(R) with alternating proximal updates, compared to C SLRA';
                    fgName      = 'fminlbfgs_prox';
                end

                use_iter = 0;
                if ~isfield(fminlbfgsData, 't_stamps'), ...
                        fminlbfgsData.t_stamps = 1:length(fminlbfgsData.f_slra_val); use_iter = 1;end


                % AGGREGATE FIGURES 
                aggrTstmps = [aggrTstmps fminlbfgsData.t_stamps(end)];
                figure(aggrMfig)
                if ~isAccSemilog
                    plot(fminlbfgsData.t_stamps, max(fminlbfgsData.M0, 0))
                else
                    logMarray = max(mean(fminlbfgsData.M0(:,:),1), ...
                        1./(abs(mean(fminlbfgsData.M0(:,:),1))+1));
                    aggrMsMin = [aggrMsMin min(logMarray)];
                    aggrMsMax = [aggrMsMax max(logMarray)]; 
                    
                    semilogy(fminlbfgsData.t_stamps, ...
                        logMarray)
                    ylim([min(aggrMsMin) 100])
                    %semilogy(fminlbfgsData.t_stamps, mean(fminlbfgsData.M0))
                end

                figure(aggrFVALfig)
                plot(fminlbfgsData.t_stamps, fminlbfgsData.f_slra_val)
            
            end

        case 4
        %% FROM PANOC 
        % VISUALIZE FVALS, SLRA, CONSTRAINTS, AND ACCURACY. ALSO COMPARE TO OPTIMAL SLRA_MEX VALUES 

            for i = 1:2

                isFMINLBFGS = i - 1;   % If using simple lbfgs or fminlbfgs

                if isFMINLBFGS, panocData   = panocFminlbfgsData; else, ...
                        panocData = panocLbfgsData;end

                use_iter    = 0;
                if ~isfield(panocData, 't_stamps'), ...
                        panocData.t_stamps = 1:length(panocData.f_slra_val); use_iter = 1;end

                % AGGREGATE FIGURES 
                aggrTstmps = [aggrTstmps panocData.t_stamps(end)];
                figure(aggrMfig)                
                if ~isAccSemilog
                    plot(panocData.t_stamps, max(panocData.M0, 0))
                else
                    logMarray = max(mean(panocData.M0(:,:),1), ...
                        1./(abs(mean(panocData.M0(:,:),1))+1));
                    aggrMsMin = [aggrMsMin min(logMarray)];
                    aggrMsMax = [aggrMsMax max(logMarray)];                     
                    semilogy(panocData.t_stamps, ...
                        logMarray)
                    ylim([min(aggrMsMin) 100])
                    % semilogy(panocData.t_stamps, mean(panocData.M0))
                end

                figure(aggrFVALfig)
                plot(panocData.t_stamps, panocData.f_slra_val)
            
            end

        case 5 
        %% FROM ALM
        plotALMdata = almData;
        almVisualize;
        fgName      = 'alm';
        set(gcf, 'Position', get(0, 'Screensize'));
        saveas(gcf,['temp\' fgName], 'png')

    end
end

% AGGREGATE FIGURES 

% ALL : 'gd_simple', 'gd_reg', 'gd_manopt', 'gd_proj', ...
      % 'fmincon_alm','fmincon_gdTrue','fmincon_gdFalse','fmincon_varpro', ...
      % 'fminlbfgs_simple', 'fminlbfgs_prox' ,'panoc_{lbfgs}', 'panoc_{fminlbfgs}')

    % AGGREAGATE FVALS
figure(aggrFVALfig)

plot(info_ident.iterinfo(1,:), info_ident.iterinfo(2,:), 'r--')
line([0 mean(aggrTstmps)], [f_slra f_slra], 'Color', 'r')

legend('gd_{reg}', 'gd_{proj}', ...
    'fmincon', ...
    'fminlbfgs_{simple}', 'fminlbfgs_{prox}', ...
    'panoc_{lbfgs}', 'panoc_{fminlbfgs}', ...
    'slra mex', 'slra optimum')

xlim([0 mean(aggrTstmps)])
xlabel('t (seconds)')

fgDetails = ['-C', num2str(statsTable.complexities(cmplx_iter)), '-m', num2str(m_in), '-p',num2str(p_out)];
fgName = ['Aggr' fgDetails '-FVALS'];


set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,['temp\' fgName], 'mfig')
saveas(gcf,['temp\' fgName], 'png')

    % AGGREAGATE ACCURACIES 
figure(aggrMfig)
line([0.001 mean(aggrTstmps)], ...
    [mean(M_slra) mean(M_slra)], 'Color', 'r')

legend('gd_{reg}', 'gd_{proj}', ...
    'fmincon', ...
    'fminlbfgs_{simple}', 'fminlbfgs_{prox}', ...
    'panoc_{lbfgs}', 'panoc_{fminlbfgs}', ...
    'slra best accuracy')

xlim([0.001 mean(aggrTstmps)])
ylim([min(aggrMsMin) 100])
xlabel('t (seconds)')
set(gca,'yscale','log')

fgName = ['Aggr' fgDetails '-Ms'];
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf,['temp\' fgName], 'mfig')        
saveas(gcf,['temp\' fgName], 'png')        
