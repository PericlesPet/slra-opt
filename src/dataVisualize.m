% Script for All Data Visualizations
% select Visual:
% 1 :   gradient descent Data
% 2 :   myFmincon        Data
% 3 :   fminlbfgs        Data 
% 4 :   panoc            Data
% 5 :   alm              Data
if isCloseAll == 1
    close all
    clc
else
    isCloseAll = 1;
end

    % Load ALL Data
all_files = dir(fullfile('data\final Data\','*.mat'));
for i = 1:numel(all_files)
    if ~exist(all_files(i).name(1:end-4)) ... 
            || strcmp(all_files(i).name(1:end-4), 'fminconData_aLM') ...
            || strcmp(all_files(i).name(1:end-4), 'almData')
        file_path = [all_files(i).folder '\' all_files(i).name];
        fprintf('Loading: %s \n',file_path);
        load(file_path)
    else 
        fprintf('Already Loaded: %s \n',all_files(i).name(1:end));
    end
end

% info_ident = 1;

for selectVisual = 1:5
    switch selectVisual

        case 1
        %% FROM GD
            gdXlim = 2;

            for selectGDesc = 1:4
                if selectGDesc == 1
                    gdescData   = gdData;    
                    gdescTitle  = 'Simple Gradient Descent';
                    fgName      = 'gd_simple';
                elseif selectGDesc == 2
                    gdescData   = gdRegData;    
                    gdescTitle  = 'Regularized Gradient Descent';
                    fgName      = 'gd_reg';
                elseif selectGDesc == 3
                    gdescData   = gdManoptData;    
                    gdescTitle  = 'ManOpt-like Gradient Descent';
                    fgName      = 'gd_manopt';
                elseif selectGDesc == 4
                    gdescData   = gdProjData;    
                    gdescTitle  = 'Stiefel-Projected Gradient Descent';
                    fgName      = 'gd_proj';
                end

                use_iter = 0;
                if ~isfield(gdescData, 't_stamps'), ...
                        gdescData.t_stamps = 1:length(gdescData.f); use_iter = 1;end
                set(0,'DefaultLegendAutoUpdate','off')

                figure
                subplot(3,1,1)
                plot(gdescData.t_stamps, gdescData.f)
                title('SLRA M(R)')
                hold on
                if ~use_iter, plot(info_ident.iterinfo(1,:), info_ident.iterinfo(2,:), 'r--'), else ...
                    plot(1:length(info_ident.iterinfo(1,:)), info_ident.iterinfo(2,:), 'r--'), end
                line([0 gdescData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')
                legend('GD', 'slra mex', 'optimum point')
                xlim([0 gdXlim])
                xlabel('t (seconds)')

                subplot(3,1,2)
                for ii = 1:size(gdescData.residual_R, 1)    
                    semilogy(gdescData.t_stamps, gdescData.residual_R(ii,:))
                    hold on
                end
                % hold on
                % semilogy(gdescData.CE(2,:))
                title('Constraints')
                if size(gdescData.residual_R, 1) == 1
                    legend('RR^T - I_N = 0 Const.')
                else
                    legend('R*Hank(p_{hat}) = 0 Const.', 'RR^T - I_N = 0 Const.')
                end
                xlim([0 gdXlim])
                xlabel('t (seconds)')

            %     subplot(2,2,3)
            %     plot(gdescData.t_stamps, gdescData.fvals)
            %     title('Fvals')

                subplot(3,1,3)
                plot(gdescData.t_stamps, max(gdescData.M, -100), 'b')
                % ylim([-200 100])
                title('Accuracy (%)')
                if size(gdescData.M, 1) == 2
                    line([0 gdescData.t_stamps(end)], ...
                        [M_slra(1) M_slra(1)], 'Color', 'r')
                    line([0 gdescData.t_stamps(end)], ...
                        [M_slra(2) M_slra(2)], 'Color', 'r')
                    legend('accuracy for y1', 'accuracy for y2', 'slra best accuracy for y1','slra best accuracy for y2')
                else
                    line([0 gdescData.t_stamps(end)], ...
                        [mean(M_slra) mean(M_slra)], 'Color', 'r')
                    legend('accuracy for fmincon', 'slra best accuracy')
                end
                xlim([0 gdXlim])
                xlabel('t (seconds)')
                ylim([-100 100])
                suptitle(gdescTitle)                
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(gcf,['temp\' fgName], 'png')
            end


        case 2
        %% FROM myFMINCON

        % %% VISUALIZE FVALS, SLRA, CONSTRAINTS, AND ACCURACY. ALSO COMPARE TO OPTIMAL SLRA_MEX VALUES
        % ----> dataVisualize 2ND PART


        for selectFminconFunc = 1:4
            if selectFminconFunc == 1
                fminconData = fminconData_aLM;
                fminconTitle = 'Matlab Fmincon for: L(x,R) s.t. RR^{T} - I = 0, compared to C SLRA';
                fgName      = 'fmincon_alm';
            elseif selectFminconFunc == 2
                fminconData = fminconData_gdTrue;
                fminconTitle = 'Matlab Fmincon for: f(x) = |p-x|^{2} s.t. both constraints, compared to C SLRA';
                fgName      = 'fmincon_gdTrue';
            elseif selectFminconFunc == 3
                fminconData = fminconData_gdFalse;
                fminconTitle = 'Matlab Fmincon for: f(x) = |p-x|^{2} s.t. both constraints, compared to C SLRA';
                fgName      = 'fmincon_gdFalse';
            elseif selectFminconFunc == 4
                fminconData = fminconData_slraVSslra;
                fminconTitle = 'Matlab Fmincon for: M(R) s.t. RR^{T} - I = 0, compared to C SLRA';
                fgName      = 'fmincon_varpro';
            end


            use_iter = 0;
            if ~isfield(fminconData, 't_stamps'), ...
                    fminconData.t_stamps = 1:length(fminconData.f_slra_val); use_iter = 1;end


            set(0,'DefaultLegendAutoUpdate','off')
            figure
            subplot(2,2,1)
            plot(fminconData.t_stamps, fminconData.f_slra_val)
            title('SLRA M(R)')
            hold on
            if ~use_iter, plot(info_ident.iterinfo(1,:), info_ident.iterinfo(2,:), 'r--'), else ...
                plot(1:length(info_ident.iterinfo(1,:)), info_ident.iterinfo(2,:), 'r--'), end
            line([0 fminconData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')
            legend('fmincon', 'slra mex', 'optimum point')
            xlabel('t (seconds)')

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
            xlabel('t (seconds)')

            subplot(2,2,3)
            plot(fminconData.t_stamps, fminconData.fvals)
            title('Fvals')
            xlabel('t (seconds)')

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
            xlabel('t (seconds)')
            ylim([-100 100])

            suptitle(fminconTitle)
            
            set(gcf, 'Position', get(0, 'Screensize'));
            saveas(gcf,['temp\' fgName], 'png')
        end

        % subplot(2,2,3)
        % plot(fminconData.f_slra_val)

        % VISUALIZE DIFFERENCE BETWEEN FVALS, SLRA, AND CONSTRAINTS (PROBABLY USELESS)
    %     
    %     for visualizeDiff = 1
    %         figure
    %         %%%%% SLRA %%%%%
    %         subplot(2,2,1)
    %         plot(fminconData.t_stamps, fminconData.f_slra_val)
    %         title('SLRA M(R)')
    %         hold on
    %         if ~use_iter, plot(info_mex.iterinfo(1,:), info_mex.iterinfo(2,:), 'r--'), else ...
    %             plot(1:length(info_mex.iterinfo(1,:)), info_mex.iterinfo(2,:), 'r--'), end
    %         line([0 fminconData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')
    %         legend('fmincon', 'slra mex', 'optimum point')
    % 
    %         %%%%% FVALS %%%%%
    %         subplot(2,2,3)
    %         plot(fminconData.t_stamps, 2*fminconData.fvals)
    %         title('Fvals')
    % 
    %         %%%%% SLRA - FVALS %%%%%
    %         subplot(2,2,2)
    %         plot(fminconData.t_stamps, fminconData.f_slra_val-2*fminconData.fvals)
    %         title('SLRA - Fvals')
    % 
    %         %%%%% CONSTRAINTS%%%%%
    %         subplot(2,2,4)
    %         for ii = 1:size(fminconData.CE, 1)    
    %             semilogy(fminconData.t_stamps, fminconData.CE(ii,:))
    %             hold on
    %         end
    %         % hold on
    %         % semilogy(fminconData.CE(2,:))
    %         title('Constraints')
    %         if size(fminconData.CE, 1) == 1
    %             legend('RR^T - I_N = 0 Const.')
    %         else
    %             legend('R*Hank(p_{hat}) = 0 Const.', 'RR^T - I_N = 0 Const.')
    %         end
    %     end



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


                set(0,'DefaultLegendAutoUpdate','off')
                figure
                subplot(2,2,1)
                plot(fminlbfgsData.t_stamps, fminlbfgsData.f_slra_val)
                title('SLRA M(R)')
                hold on
                if ~use_iter, plot(info_ident.iterinfo(1,:), info_ident.iterinfo(2,:), 'r--'), else ...
                    plot(1:length(info_ident.iterinfo(1,:)), info_ident.iterinfo(2,:), 'r--'), end
                line([0 fminlbfgsData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')
                legend('fmincon', 'slra mex', 'optimal slra')
                ylabel('M(R)')
                xlabel('t (seconds)')

                subplot(2,2,2)
                for ii = 1:size(fminlbfgsData.CE, 1)    
                    semilogy(fminlbfgsData.t_stamps, fminlbfgsData.CE(ii,:))
                    hold on
                end
                % hold on
                % semilogy(lbfgsData.CE(2,:))
                title('Constraints')
                if size(fminlbfgsData.CE, 1) == 1
                    legend('RR^T - I_N = 0 Const.')
                else
                    legend('R*Hank(p_{hat}) = 0 Const.', 'RR^T - I_N = 0 Const.')
                end
                xlabel('t (seconds)')

                subplot(2,2,3)
                plot(fminlbfgsData.t_stamps, fminlbfgsData.fvals)
                title('Fvals')
                ylabel('f(x)')
                xlabel('t (seconds)')

                subplot(2,2,4)
                plot(fminlbfgsData.t_stamps, max(fminlbfgsData.M0, -100), 'b')
                % ylim([-200 100])
                title('Accuracy (%)')
                if size(fminlbfgsData.M0, 1) == 2
                    line([0 fminlbfgsData.t_stamps(end)], ...
                        [M_slra(1) M_slra(1)], 'Color', 'r')
                    line([0 fminlbfgsData.t_stamps(end)], ...
                        [M_slra(2) M_slra(2)], 'Color', 'r')
                    legend('accuracy for y1', 'accuracy for y2', 'slra best accuracy for y1','slra best accuracy for y2')
                else
                    line([0 fminlbfgsData.t_stamps(end)], ...
                        [mean(M_slra) mean(M_slra)], 'Color', 'r')
                    legend('accuracy for fmincon', 'slra best accuracy')
                end
                ylabel('accuracy (%)')
                xlabel('t (seconds)')
                ylim([-100 100])

                suptitle(fminlbfgsTitle)
                
                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(gcf,['temp\' fgName], 'png')

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


                set(0,'DefaultLegendAutoUpdate','off')
                figure
                subplot(2,2,1)
                plot(panocData.t_stamps, panocData.f_slra_val)
                title('SLRA M(R)')
                hold on
                if ~use_iter, plot(info_ident.iterinfo(1,:), info_ident.iterinfo(2,:), 'r--'), else ...
                    plot(1:length(info_ident.iterinfo(1,:)), info_ident.iterinfo(2,:), 'r--'), end
                line([0 panocData.t_stamps(end)], [f_slra f_slra], 'Color', 'r')

                if ~isFMINLBFGS, legend('panoc_{lbfgs}', 'slra mex', 'optimal slra'), else ...
                        legend('panoc_{fminlbfgs}', 'slra mex', 'optimal slra'),end


                ylabel('M(R)')
                xlabel('t (seconds)')

                subplot(2,2,2)
                for ii = 1:size(panocData.CE, 1)    
                    semilogy(panocData.t_stamps, panocData.CE(ii,:))
                    hold on
                end
                % hold on
                % semilogy(lbfgsData.CE(2,:))
                title('Constraints')
                if size(panocData.CE, 1) == 1
                    legend('RR^T - I_N = 0 Const.')
                else
                    legend('R*Hank(p_{hat}) = 0 Const.', 'RR^T - I_N = 0 Const.')
                end
                xlabel('t (seconds)')

                subplot(2,2,3)
                plot(panocData.t_stamps, panocData.fvals)
                title('Fvals')
                ylabel('f(x)')
                xlabel('t (seconds)')

                subplot(2,2,4)
                plot(panocData.t_stamps, max(panocData.M0, -100), 'b')
                % ylim([-200 100])
                title('Accuracy (%)')
                if size(panocData.M0, 1) == 2
                    line([0 panocData.t_stamps(end)], ...
                        [M_slra(1) M_slra(1)], 'Color', 'r')
                    line([0 panocData.t_stamps(end)], ...
                        [M_slra(2) M_slra(2)], 'Color', 'r')
                    legend('accuracy for y1', 'accuracy for y2', 'slra best accuracy for y1','slra best accuracy for y2')
                else
                    line([0 panocData.t_stamps(end)], ...
                        [mean(M_slra) mean(M_slra)], 'Color', 'r')
                    legend('accuracy for panoc_{lbfgs}', 'slra best accuracy')
                end
                ylabel('accuracy (%)')
                xlabel('t (seconds)')
                ylim([-100 100])

                if isFMINLBFGS 
                    suptitle('PANOC with FMINLBFGS, compared to C SLRA')
                    fgName      = 'panoc_fminlbfgs';
                else
                    suptitle('PANOC with simple LBFGS, compared to C SLRA')
                    fgName      = 'panoc_lbfgs';
                end

                set(gcf, 'Position', get(0, 'Screensize'));
                saveas(gcf,['temp\' fgName], 'png')


                % subplot(2,2,3)
                % plot(lbfgsData.f_slra_val)
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









        