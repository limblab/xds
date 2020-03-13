function [trial_spike_counts, trial_EMG, trial_force, trial_kin] = get_rewarded_trials(xds)
suc_timetable = get_trial_time_table(xds, 'R');
j = 1;
for i = 1:length(suc_timetable)
    if suc_timetable(i,2) - suc_timetable(i,1)<0.5
        continue;
    end
    if isnan(suc_timetable(i,1)) || isnan(suc_timetable(i,2))
       continue;
    end
    temp = find((xds.time_frame >= suc_timetable(i,1))&(xds.time_frame <= suc_timetable(i,2)));
    trial_spike_counts{j,1} = xds.spike_counts(temp, :);
    if xds.has_EMG == true
       trial_EMG{j,1} = xds.EMG(temp, :);
    else
       trial_EMG = 0;
    end
    if xds.has_force == true
       trial_force{j,1} = xds.force(temp, :);
    else
       trial_force = 0;
    end
    if xds.has_kin == true
       trial_kin{j,1} = xds.kin_p(temp, :);
       trial_kin{j,2} = xds.kin_v(temp, :);
       trial_kin{j,3} = xds.kin_a(temp, :);
    else
       trial_kin = 0;
    end
    j = j+1;
end

end

