clc
clear
params = struct( ...
    'monkey_name','Pop', ...
    'array_name','M1', ...
    'task_name','FR', ...
    'ran_by','MP', ...
    'lab',1, ...
    'bin_width',1/30, ...
    'sorted',1, ...
    'requires_raw_emg',1);

base_dir = 'D:\Kevin\';
file_name = 'Pop_20201103_FR_001.nev';
map_dir = 'Z:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\Implant_2020_01\6250-002086\';
map_name = 'SN 6250-002086.cmp';
save_dir = 'D:\Kevin\';
xds = raw_to_xds(base_dir, file_name, map_dir, map_name, params);
save_name = strcat(save_dir, file_name);
save(strcat(save_name, '.mat'),'xds');