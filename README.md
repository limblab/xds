# xds (cross-platform data structure), as a wrap for `cds` to be used in Python

The latest updates are in the branch "wireless_EMG", which also supports EMGs collected using wired systems. Will merge with the main branch soon.

## Overview
`cds` (Classy data structure) is the standard data structure to load data collected in lab. Currently, Python packages including scipy and h5py are not friendly to MATLAB classes and tables, so here we create this wrap to break those high-level MATLAB structures into simpler components. 

## How to create xds files from .nev and .nsx files?

* Add all necessary folders to your Matlab path:

	- ClassyDataAnalysis (https://github.com/limblab/ClassyDataAnalysis) and all subdirectories (don't need the .git directories, but feel free to remove them)
	- put xds_matlab in your MATLAB path

* You can write a script as below to load your data, or simply open the file xds_matlab/example_load_single_file.m, and change the params and file destinations.


```
clc
clear
params = struct( ...
    'monkey_name','Pop', ...
    'array_name','M1', ...
    'task_name','WM', ...
    'ran_by','KLB', ...
    'lab',1, ...
    'bin_width',0.001, ...
    'sorted',0, ...
    'requires_raw_emg',1);

base_dir = 'Z:\data\Pop_18E3\CerebusData\20200311\';
file_name = 'Pop_20200311_WM_CO_001';
map_dir = 'Z:\limblab\lab_folder\Animal-Miscellany\Pop_18E3\Array Map Files\Implant_2020_01\6250-002086\';
map_name = 'SN 6250-002086.cmp';
save_dir = '.\';
xds = raw_to_xds(base_dir, file_name, map_dir, map_name, params);
save_name = strcat(save_dir, file_name);
save(strcat(save_name, '.mat'),'xds');
```

## How to create xds files from many .nev and .nsx files?

* Sometimes you may need to convert a series of files to xds format togother. Here, you can start with putting all `.nev` and `.nsx` in one folder.

* Then, you can create a script like this:
 ```
clc
clear
params = struct( ...
    'monkey_name','Jango', ...
    'array_name','M1', ...
    'task_name','WF', ...
    'ran_by','SN', ...
    'lab',1, ...
    'bin_width',0.001,...
    'sorted',0,...
    'requires_raw_emg',1);

base_dir = 'Z:\limblab\lab_folder\Projects\darpa\DS18(Jango_2015)\nev\';
map_dir = 'Z:\limblab\lab_folder\Projects\darpa\array_map_files\Jango_right_M1\';
map_name = 'SN6250-000945.cmp';
save_dir = '.\'; 
open_file = strcat(base_dir, '*.nev');
file = dir(open_file);
for ii = 1:length(file)
    file_name_in_list = file(ii).name(1:end-4);
    disp(file_name_in_list);
    xds = raw_to_xds(base_dir, file_name_in_list, map_dir, map_name, params);
    save_file = strcat(file_name_in_list, '.mat');
    save(strcat(save_dir, save_file), 'xds');
    clear xds
end
```
* Finally, you will find a series of files in your `save_dir`

## How to use xds in Python?
Check readme.md in folder `xds_python`.
