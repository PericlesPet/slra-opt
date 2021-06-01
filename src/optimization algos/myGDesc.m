function [logdata, dataOptId, f_log, minf_log] = myGDesc(gdInput, opt, R_slramex)
    %% Parameters
    clear logdata;
    Rin     = gdInput.Rin;
    maxIter = gdInput.maxIter;
    gamma   = gdInput.gamma;
    reg     = gdInput.reg;
    obj     = gdInput.obj;
    Ropt    = gdInput.Ropt;
    extCond = gdInput.extCond;
    sysAcc  = gdInput.sysAccuracy;
    % last_x_elements = 10000;
    % last_x_elements = maxIter;
    % maxIter = ceil(maxIter / last_x_elements) * last_x_elements;
    f_log = zeros(maxIter, 1);
    minf_log = zeros(ceil(maxIter/100),1);
    min_f = inf;
    % Regularizer Parameter
    if reg
        mu = opt.g;
    end

    modVal = 25;
    % Gradient Update Parameters
    constant_thing = 0.1;   % To avoid division by zero
    stepsize_factor = 0.01;  % Initial value of stepsize parameter
    stepsize_update = 20;  % Reduce stepsize every X steps
    % If the optimal value changes less than (tolPcntChange)% 
    % across 100 iterations, stop the descent
    tolPcntChange = 1;      

    %%
    tic
    for i = 1:maxIter

        %% Set F
        f = slra_mex_obj('func', obj, Rin);
        % Rin;
        if reg
            regularizer = norm(Rin * Rin' - eye(size(Rin*Rin')),'fro')^2;
            f_reg = f + mu * regularizer;
        end

        if i == 2
            if reg
                min_f = f_reg;
            else
                min_f = f ; 
            end
        end


        %% Set G
        if exist('g') prev_dir = sign(g); else prev_dir = zeros(size(Rin)); end
        g = slra_mex_obj('grad', obj, Rin);
        curr_dir = sign(g);
        if reg
            regularizer_grad =  2 *(Rin * Rin' - eye(size(Rin*Rin')))*Rin;
            g_reg = g + mu * regularizer_grad ;
        end


        %% Every modVal iterations ...
        printData = 1;
        if (mod(i-1,modVal) == 0) && (printData == 1)
            %% Data Logging
            if exist('logdata'), t_stamp = logdata.t_stamps(end) + toc; ...
            else, t_stamp = toc; end
            residual_R = norm(Rin*Rin'-eye(size(Rin,1)));
            %     dir_diff    = sign(R-Rin) - sign(-g);
            %         id = mod(i-1, last_x_elements)+1;
            id = (i-1) / modVal + 1;
            f_log(i) = f; 
            logdata.t_stamps(id)      = t_stamp;
            logdata.f(id)             = f;
            logdata.g(:,id)           = g(:);
            logdata.residual_R(id)    = residual_R;
            logdata.Rin(:,id)         = Rin(:);
            if reg
                logdata.f_reg(id)         = f_reg;
                logdata.g_reg(:,id)       = g_reg(:);
                %         logdata(id).pcntChange  = norm(gamma*g_reg)/norm(Rin);
                %         logdata(id).normG_reg   = norm(gamma*g_reg);
            else
                %         logdata(id).pcntChange  = norm(gamma*g)/norm(Rin);
                %         logdata(id).normG       = norm(gamma*g);
            end
            logdata.normRin(id)          = norm(Rin);
            [~, logdata.M(:,id)]         = sysAcc(Rin);
            %     logdata(id).R_opt       = R;
            %     logdata(id).R_opt_dist  = R - Rin;
            %     logdata(id).dir_diff    = dir_diff;
            %     logdata(id).R_cond      = cond(Rin);
            %     logdata(id).gdDir_diff  = prev_dir - curr_dir;


            %% Log minimum
            if reg
                if f_reg <= min_f
                    dataOptId = id;
                    min_f = f_reg;  
                end
            else
                if f<=min_f
                    dataOptId = id;
                    min_f = f;
                end
            end

            %% Condition + Print
    %         logdata
            fprintf('ITERATION: %d, minF = %f, id = %d, t_stamp = %f \n', i, min_f, id, logdata.t_stamps(id));
            minf_log(id) = min_f;
            if (id >= 2 && extCond)
                %         abs((f_log(i-1)/f_log(i - 101)-1)*100)
                abs((minf_log(i/modVal)/minf_log(i/modVal-1)-1)*100)        
                if abs((minf_log(i/modVal)/minf_log(i/modVal-1)-1)*100) < tolPcntChange
                    break
                end
            end
            tic
        else
            dataOptId = 1;
        end

        %% Update R
        stepsize_param = stepsize_factor/ceil(i/stepsize_update);

        if reg    
            Rin = Rin - stepsize_param * (gamma * g_reg) / (norm(gamma*g_reg)+constant_thing) ;
        else
            Rin = Rin - stepsize_param * gamma * g / (norm(gamma*g) + constant_thing)  ;
        end

    end

    logdata.gamma   = gamma;
    logdata.reg     = reg;
    logdata.updateParams = [constant_thing;stepsize_factor;stepsize_update];
    logdata.comment = 'Update Params: constant_thing, stepsize_factor, stepsize_update';


