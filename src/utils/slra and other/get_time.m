function [finishInfo] = get_time(data)

if isfield(data, 'fvals')
    fvals = data.fvals;
    Ms    = data.M0;
else
    fvals = data.f;
    Ms    = data.M;
end
tstamps = data.t_stamps;

ma_points = 6;
pcnt_target = 0.01;


pcntgs = ((fvals(1:end-1) - fvals(2:end))./fvals(1:end-1))'*100;
clear mavg

for i = 1:length(pcntgs)
    % MA index
    ma_i = mod(i-1,ma_points)+1;
    % Moving Average Array
    mavg(ma_i) = pcntgs(i);
    % Mean of the Array
    avg_pcnt = mean(mavg);
    
    finish_index = i+1;
    t_finish = tstamps(finish_index);
    acc_finish = Ms(:,finish_index);

    if avg_pcnt <= pcnt_target
        break;
    end

end

finishInfo.t     = t_finish;
finishInfo.M     = acc_finish;
finishInfo.Mavg  = mean(acc_finish);
finishInfo.index = finish_index;

end