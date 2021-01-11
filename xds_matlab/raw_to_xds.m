function xds = raw_to_xds(file_dir, file_name, map_dir, map_name, params)
data_file = strcat(file_dir, file_name);
map_file = strcat(map_dir, map_name);
monkey_name = params.monkey_name;
array_name = params.array_name;
task_name = params.task_name;
lab = params.lab;
ran_by = params.ran_by;
sorted = params.sorted;

% Here are some fields for raw data
requires_raw_emg = 0;
requires_raw_force = 0;
if isfield(params,'requires_raw_emg')
    requires_raw_emg = params.requires_raw_emg;
end
if isfield(params,'requires_raw_force')
    requires_raw_force = params.requires_raw_force;
end

cds=commonDataStructure();
cds.file2cds(data_file,['array', array_name],...
            ['monkey', monkey_name],lab,'ignoreJumps',['task', task_name], ...
            ['ranBy', ran_by], ['mapFile', map_file]);

bin_width = params.bin_width;
ex = experiment; 
ex.meta.hasEmg = cds.meta.hasEmg; 
ex.meta.hasUnits = true;
ex.meta.hasTrials = true; 
ex.meta.hasForce = cds.meta.hasForce;
ex.meta.hasKinematics = cds.meta.hasKinematics;
ex.meta.hasLfp = cds.meta.hasLfp;
if isempty(cds.kin) && ~isempty(cds.cursor)
    ex.meta.hasKinematics = 1;
end
    
ex.binConfig.filterConfig.sampleRate = 1/bin_width;
ex.binConfig.filterConfig.poles = 2;
ex.binConfig.filterConfig.cutoff = 10;

ex.firingRateConfig.sampleRate = 1/bin_width;
ex.firingRateConfig.method = 'bin';
ex.addSession(cds);
ex.calcFiringRate;

ex.binConfig.include(1).field = 'units';
field_ind = 2;
if ex.meta.hasEmg == true
   ex.emg.processDefault;
   ex.binConfig.include(field_ind).field = 'emg';
   field_ind = field_ind + 1;
end

if ex.meta.hasKinematics == true
   ex.binConfig.include(field_ind).field = 'kin';
   field_ind = field_ind + 1;
end

if ex.meta.hasForce == true
   ex.binConfig.include(field_ind).field = 'force';
end

% if ex.meta.hasLfp == true
%    ex.binConfig.include(field_ind).field = 'lfp';
% end

%^The above block of text is sometimes necessary if you're converting files
%from monkeys that have multiple implants. Such a file can have only
%units and LFP with no EMG/kinematics/cursor data, because the
%EMG/kinematics/cursor data is in the file recorded from the other array;
%trying to create an XDS with only units causes the experiment processing to
%fail, and so that requires the inclusion of LFP in the experiment processing. 
ex.binData;

xds.meta = cds.meta;
xds.bin_width = bin_width;
xds.time_frame = ex.bin.data.t;
xds.has_EMG = xds.meta.hasEmg;
xds.has_force = xds.meta.hasForce;
xds.has_kin = xds.meta.hasKinematics;
xds.sorted = sorted;

% units
elec_mask = zeros(length(cds.units), 1);
for i = 1:length(cds.units)
    if strfind(cds.units(i).label,'elec') == 1
        elec_mask(i) = 1;
    end
end

if sorted == 0
    good_id = find(elec_mask == 1);
    for i = 1:length(good_id)
        xds.unit_names{1,i} = cds.units(good_id(i)).label;
        xds.spikes{1,i} = cds.units(good_id(i)).spikes.ts;
        %xds.spike_waveforms{1,i} = cds.units(good_id(i)).spikes.wave;
    end
    [~,binnedUnitMask] = ex.bin.getUnitNames;
    bad_id = find(elec_mask == 0);
    temp1 = find(binnedUnitMask == 1);
    for i = 1:length(bad_id)
        binnedUnitMask(temp1(bad_id(i))) = 0;
    end
elseif sorted == 1
    disp('Sorted version will be served later');
end

xds.spike_counts = ex.bin.data{:,binnedUnitMask}*bin_width;
if ex.meta.hasEmg == true
   emgMask = ~cellfun(@(x)isempty(strfind(x,'EMG')),ex.bin.data.Properties.VariableNames);
   emgNames = ex.bin.data.Properties.VariableNames(emgMask);
   xds.EMG = ex.bin.data{:,emgMask};
   xds.EMG_names = emgNames;
   if requires_raw_emg
       raw_EMG_table = table2array(cds.emg);
       xds.raw_EMG = raw_EMG_table(:,2:end);
       xds.raw_EMG_time_frame = raw_EMG_table(:,1);
   else
       xds.raw_EMG = [];
       xds.raw_EMG_time_frame = [];
   end
end
if ex.meta.hasForce == true
   fxMask = ~cellfun(@(x)isempty(strfind(x,'fx')),ex.bin.data.Properties.VariableNames);
   fyMask = ~cellfun(@(x)isempty(strfind(x,'fy')),ex.bin.data.Properties.VariableNames);
   xds.force(:, 1) = ex.bin.data{:,fxMask};
   xds.force(:, 2) = ex.bin.data{:,fyMask};
   if requires_raw_force
       raw_force_table = table2array(cds.force);
       xds.raw_force = raw_force_table(:,2:end);
       xds.raw_force_time_frame = raw_force_table(:,1);
   else
       xds.raw_force = [];
       xds.raw_force_time_frame = [];
   end
end
if ex.meta.hasKinematics == true
   xMask = ~cellfun(@(x)isempty(strfind(x,'x')),ex.bin.data.Properties.VariableNames);
   yMask = ~cellfun(@(x)isempty(strfind(x,'y')),ex.bin.data.Properties.VariableNames);
   temp = find(xMask==1);
   for i = 2:length(temp)
       xMask(temp(i)) = 0; 
   end
   temp = find(yMask==1);
   for i = 2:length(temp)
       yMask(temp(i)) = 0; 
   end
   vxMask = ~cellfun(@(x)isempty(strfind(x,'vx')),ex.bin.data.Properties.VariableNames);
   vyMask = ~cellfun(@(x)isempty(strfind(x,'vy')),ex.bin.data.Properties.VariableNames);
   axMask = ~cellfun(@(x)isempty(strfind(x,'ax')),ex.bin.data.Properties.VariableNames);
   ayMask = ~cellfun(@(x)isempty(strfind(x,'ay')),ex.bin.data.Properties.VariableNames);
   xds.kin_p(:, 1) = ex.bin.data{:, xMask};
   xds.kin_p(:, 2) = ex.bin.data{:, yMask};
   xds.kin_v(:, 1) = ex.bin.data{:, vxMask};
   xds.kin_v(:, 2) = ex.bin.data{:, vyMask};
   xds.kin_a(:, 1) = ex.bin.data{:, axMask};
   xds.kin_a(:, 2) = ex.bin.data{:, ayMask};
   if isempty(cds.kin) && ~isempty(cds.cursor)
       xds.cursor_p = xds.kin_p;
       xds.cursor_v = xds.kin_v;
       xds.cursor_a = xds.kin_a;
       xds = rmfield(xds,'kin_p');
       xds = rmfield(xds,'kin_v');
       xds = rmfield(xds,'kin_a');
       xds.has_kin = 0;
       xds.has_cursor = 1;
   else
       xds.has_cursor = 0;
   end
end   

% trial information
xds.trial_info_table_header = fieldnames(cds.trials);
xds.trial_info_table = table2cell(cds.trials);

xds.trial_gocue_time = deal_trial_info('goCue', cds);
xds.trial_start_time = deal_trial_info('startTime', cds);
xds.trial_end_time = deal_trial_info('endTime', cds);
xds.trial_result = deal_trial_info('result', cds);
xds.trial_target_dir = deal_trial_info('tgtDir', cds);
xds.trial_target_corners = deal_trial_info('Corners', cds);

clear cds
clear ex
end

function trial_info = deal_trial_info(str,cds)
trial_mask = ~cellfun(@(x)isempty(strfind(x,str)),cds.trials.Properties.VariableNames);
if sum(trial_mask) == 0
    disp('Something is wrong with the trial table');
    trial_info = 0;
else 
    trial_info = cds.trials{:, trial_mask};
end
end