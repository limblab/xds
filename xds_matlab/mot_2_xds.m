
function [joint_angle_time_frame, joint_angles, joint_names] = mot_2_xds(video_sync, file_dir, file_name)

%% Load the .mot file

xds_mot = readtable(strcat(file_dir, file_name, '.mot'), 'Filetype', 'text');
disp('Loading mot file:')

mot_var_names = xds_mot.Properties.VariableNames';

% Load the joint angle time frames
time_frame_idx = find(contains(mot_var_names, 'time'));
joint_angle_time_frame = xds_mot.(time_frame_idx);

joint_names = struct([]);
joint_angless = struct([]);
cc = 1;
for jj = 1:length(mot_var_names)
    joint_title = mot_var_names{jj,1};
    joint_names{cc} = '';
    if contains(joint_title, 'Thorax') || contains(joint_title, 'clavicular')
        continue
    end
    if contains(joint_title, 'unro') || contains(joint_title, 'shoulder') || contains(joint_title, 'prox')
        continue
    end
    % Which finger
    if contains(joint_title, '1')
        joint_names{cc} = 'Thumb_';
    elseif contains(joint_title, '2')
        joint_names{cc} = 'Index_';
    elseif contains(joint_title, '3')
        joint_names{cc} = 'Middle_';
    elseif contains(joint_title, '4')
        joint_names{cc} = 'Ring_';
    elseif contains(joint_title, '5')
        joint_names{cc} = 'Pinky_';
    elseif contains(joint_title, 'wr')
        joint_names{cc} = 'Wrist_';
    end
    if isempty(joint_names{cc})
        continue
    end
    % Which joint
    if contains(joint_title, 'dip')
        joint_names{cc} = strcat(joint_names{cc}, 'DIP_');
    elseif contains(joint_title, 'pip')
        joint_names{cc} = strcat(joint_names{cc}, 'PIP_');
    elseif contains(joint_title, 'mcp')
        joint_names{cc} = strcat(joint_names{cc}, 'MCP_');
    elseif contains(joint_title, '_ip')
        joint_names{cc} = strcat(joint_names{cc}, 'IP_');
    elseif contains(joint_title, 'cmc')
        joint_names{cc} = strcat(joint_names{cc}, 'CMC_');
    end
    % Which angle
    if contains(joint_title, 'e_f') || contains(joint_title, 'f_e')
        joint_names{cc} = strcat(joint_names{cc}, 'Flex_Ext');
    elseif contains(joint_title, 'ad_ab')
        joint_names{cc} = strcat(joint_names{cc}, 'Add_Abd');
    elseif contains(joint_title, 'rd_ud')
        joint_names{cc} = strcat(joint_names{cc}, 'RadDev_UlnDev');
    elseif contains(joint_title, 'sup_pro')
        joint_names{cc} = strcat(joint_names{cc}, 'Sup_Pro');
    elseif contains(joint_title, 'opp')
        joint_names{cc} = strcat(joint_names{cc}, 'Opp');
    end

    joint_angless{cc} = xds_mot.(jj);
    cc = cc + 1;

end

joint_angles = zeros(length(joint_angle_time_frame), length(joint_angless));
for jj = 1:length(joint_angless)
    joint_angles(:,jj) = joint_angless{jj};
end

%% Correct joint angle time frames
% Initialization offset
offset_t = video_sync.t(find(video_sync.video_sync > 3000, 1));
joint_angle_time_frame = joint_angle_time_frame + offset_t;











