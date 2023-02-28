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

Monkey_Hand = 'Left';
TgtHold = 0.5;

params = struct( ...
    'monkey_name', 'Pop', ...
    'array_name', 'M1', ...
    'task_name', 'multi_gadget', ... % WS, multi_gadget, FR, WB, etc.
    'ran_by', 'HP', ...
    'lab', 1, ...
    'bin_width', 0.001,...
    'sorted', 1,...
    'requires_raw_emg', 1,...
    'save_waveforms', 1);

file_dir = 'C:\Users\rhpow\Documents\Work\Northwestern\Monkey_Data\Pop\20210617\';
map_dir = 'C:\Users\rhpow\Documents\Work\Northwestern\Monkey_Data\Pop\';
map_name = 'SN 6250-002339';
save_dir = 'C:\Users\rhpow\Documents\Work\Northwestern\Monkey_Data\Pop\20210617\';
open_file = strcat(file_dir, '*.ccf');
file = dir(open_file);

for ii = 1:length(file)
    file_name = file(ii).name(1:end-4);
    disp(file_name);
    xds = raw_to_xds(file_dir, file_name, map_dir, map_name, params);
    %[xds] = CalculateNonLinearEnergy(xds);
    xds.meta.hand = Monkey_Hand;
    xds.meta.TgtHold = TgtHold;
    if params.sorted == 1
        save_file = strcat(file_name, '-s', '.mat');
    else
        save_file = strcat(file_name, '.mat');
    end
    save(strcat(save_dir, save_file), 'xds', '-v7.3');
    clear xds
end
```
* Finally, you will find a series of files in your `save_dir`

## How to use xds in Python?
Check readme.md in folder `xds_python`.
