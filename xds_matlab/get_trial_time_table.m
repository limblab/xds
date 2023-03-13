function [trial_time_table] = get_trial_time_table(xds, trial_type, start_event, time_before_start, end_event, time_after_end)

% For trial_type, 'R': rewarded trials, 'A': aborted trials, 'F': failed
temp = find(xds.trial_result == trial_type);
trial_time_table = zeros(length(temp),2);

for i =1:size(temp,1)
    idx = temp(i);
    if start_event == "start_time"
        time_start = xds.trial_start_time(idx);
    elseif start_event == "gocue_time"
        time_start = xds.trial_gocue_time(idx);
    elseif start_event == "force_onset_time"
        time_start = xds.trial_force_onset_time(idx);
    elseif start_event == "end_time"
        time_start = xds.trial_end_time(idx);
    end

    if end_event == "start_time"
        time_end = xds.trial_start_time(idx);
    elseif end_event == "gocue_time"
        time_end = xds.trial_gocue_time(idx);
    elseif end_event == "force_onset_time"
        time_end = xds.trial_force_onset_time(idx);
    elseif end_event == "end_time"
        time_end = xds.trial_end_time(idx);
    end

    t1 = time_start - time_before_start;
    t2 = time_end + time_after_end;

    if t2 > t1
        trial_time_table(i,:) = [t1 t2];
    else
        trial_time_table(i,:) = [];
%         print('The timing with trial No. %d is not right.'%(n))
    end
end
end

