clc
clear
params = struct( ...
    'monkey_name','Jango', ...
    'array_name','M1', ...
    'task_name','WS', ...
    'ran_by','XM', ...
    'lab',1, ...
    'bin_width',0.001,...
    'sorted',0);

base_dir = 'F:\data\';
map_dir = 'Z:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\Implant_2020_01\6250-002086\';
map_name = 'SN 6250-002086.cmp';
save_dir = 'F:\'; 
open_file = strcat(base_dir, '*.nev');
file = dir(open_file);
for i = 2:2
    file_name_in_list = file(i).name(1:end-4);
    disp(file_name_in_list);
    xds = raw_to_xds(base_dir, file_name_in_list, map_dir, map_name, params);
    save_file = strcat(file_name_in_list, '.mat');
    save(strcat(save_dir, save_file), 'xds');
    clear xds
end

    
    
    
    
    
    
    



