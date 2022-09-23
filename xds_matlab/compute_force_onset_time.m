function idx_onset = compute_force_onset_time(trial_force)

ch = 1;
thr = 0.4;

idx_onset = zeros(length(trial_force), 1);
for i = 1:length(trial_force)
    each = trial_force{i};
%     f = each(:, ch);
    f = sqrt(each(:, 1).^2 + each(:, 1).^2);
    df = diff(f);
    temp = find(df >= thr*max(df));
    
    % Get first index after 250 ms 
    temp_idx = find(temp>6);
    
    if isempty(temp) == 0
%         idx_onset(i) = temp(1);
        idx_onset(i) = temp(temp_idx(1));
    end
end


end

