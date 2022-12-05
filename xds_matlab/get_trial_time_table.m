function [trial_time_table] = get_trial_time_table(xds, trial_type, start_time)
% For trial_type, 'R': rewarded trials, 'A': aborted trials, 'F': failed
% trials
temp = find(xds.trial_result == trial_type);
trial_time_table = zeros(length(temp),2);
for i =1:size(temp,1)
    idx = temp(i);
    if strcmp(start_time,'start_time')
        trial_time_table(i,:) = [xds.trial_start_time(idx) xds.trial_end_time(idx)];
    elseif strcmp(start_time, 'gocue_time')
        trial_time_table(i,:) = [xds.trial_gocue_time(idx) xds.trial_end_time(idx)];
    end
end
end

