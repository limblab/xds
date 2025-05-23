import numpy as np
import math
from xds import lab_data
import scipy.stats as stats
import pickle
import matplotlib.pyplot as plt
from matplotlib import rcParams
import seaborn as sns

def sort_trials_with_bump_dir(xds, trial_options, trial_type, time_params):
    """
    xds: the struct
    trial_options: 'both', 'bump', 'no bump'
    trial_type: 'R' or 'F'
    """
    start_event = time_params['start_event']
    ahead = time_params['time_before_start']
    end_event = time_params['end_event']
    end_time_offset = time_params['time_after_end']
    #---- print reach direction list ----#
    tgt_dir_list = list(set(xds.trial_target_dir))
    print('Reaching direction includes: %s'%(tgt_dir_list))
    #---- get the indices for rewarded or failed trials ----#
    trial_type_idx = np.where(xds.trial_result == trial_type)[0]
    #---- get spike count data and cursor data ----#
    spike_counts = xds.get_trials_data_spike_counts(trial_type, start_event, ahead, end_event, end_time_offset)
    curs_p = xds.get_trials_data_kin(trial_type, start_event, ahead, end_event, end_time_offset)[0]
    
    #---- sort trials according to direction and with or without bump ----#
    trial_spike_counts, trial_curs_p = [], []
    if trial_options == 'both':
        for direction in tgt_dir_list:
            target_dir_idx = np.where(xds.trial_target_dir == direction)
            b = trial_type_idx[np.isin(trial_type_idx, target_dir_idx)]
            c = sorted([np.where(trial_type_idx == each)[0][0] for each in b])
            trial_spike_counts.append([spike_counts[i] for i in c])
            trial_curs_p.append([curs_p[i] for i in c])
    elif trial_options == 'bump':
        #---- get indices for trials with a bump ----#
        bump_idx = np.where(np.isnan(xds.bump_time)!=True)[0]
        for direction in tgt_dir_list:
            target_dir_idx = np.where(xds.trial_target_dir == direction)[0]
            b = set(trial_type_idx)&set(target_dir_idx)&set(bump_idx)
            c = sorted([np.where(trial_type_idx == each)[0][0] for each in b])
            trial_spike_counts.append([spike_counts[i] for i in c])
            trial_curs_p.append([curs_p[i] for i in c])
    elif trial_options == 'no bump':
        #---- get indices for trials without a bump ----#
        bump_idx = np.where(np.isnan(xds.bump_time)==True)[0]
        for direction in tgt_dir_list:
            target_dir_idx = np.where(xds.trial_target_dir == direction)[0]
            b = set(trial_type_idx)&set(target_dir_idx)&set(bump_idx)
            c = sorted([np.where(trial_type_idx == each)[0][0] for each in b])
            trial_spike_counts.append([spike_counts[i] for i in c])
            trial_curs_p.append([curs_p[i] for i in c])
    #---- return the data ----#
    trial_curs_p_, trial_spike_counts_ = [], []
    temp = [len(each) for each in sum(trial_curs_p, [])]
    LEN = stats.mode(temp)[0]-1
    for condition in zip(trial_curs_p, trial_spike_counts):
        temp = [each[:LEN, :] for each in condition[0]]
        trial_curs_p_.append(temp)
        temp = [each[:LEN, :] for each in condition[1]]
        trial_spike_counts_.append(temp)
    return trial_spike_counts_, trial_curs_p_

def vaf(x,xhat):
    x = x - x.mean(axis=0)
    xhat = xhat - xhat.mean(axis=0)
    return (1-(np.sum(np.square(x - xhat))/np.sum(np.square(x))))

def flatten_list(X):
    """
    Converting list containing multiple ndarrays into a large ndarray
    X: a list
    return: a numpy ndarray
    """
    n_col = np.size(X[0],1)
    Y = np.empty((0, n_col))
    for each in X:
        Y = np.vstack((Y, each))
    return Y

def flatten_list_3d(X):
    n_c1 = np.size(X[0], 1)
    n_c2 = np.size(X[0], 2)
    Y = np.empty((0, n_c1, n_c2))
    for each in X:
        Y = np.vstack((Y, each))
    return Y

def flatten_list_list(my_list):
    for i in range(len(my_list)):
        my_list[i] = flatten_list(my_list[i])
    my_list = flatten_list(my_list)
    return my_list  

def xds_get_elec_idx(dataset, elec_num):
    """
    To get the idx of electrodes specified by elec_num
    dataset: xds structure
    elec_num: a list containing the number of bad channels
    """
    idx = []
    for each in elec_num:
        if 'elec'+str(each) in dataset.unit_names:
            temp = dataset.unit_names.index('elec'+str(each))
            idx.append(temp)
    return idx

def xds_del_bad_chs(dataset, elec_num):
    """
    To get rid of everything about the bad channels from my_cage_data
    """
    idx = xds_get_elec_idx(dataset, elec_num)
    for idx in sorted(idx, reverse=True):
        del(dataset.unit_names[idx])
        del(dataset.spikes[idx])
    dataset.spike_counts = np.delete(dataset.spike_counts, idx, axis = 1)
    return dataset

def get_sort_seq(spikes_list, offset):
    idx = [np.where(each>offset)[0][0] for each in spikes_list]
    t = np.flip(np.argsort([spikes_list[i][each] for i, each in enumerate(idx)]))
    return t
    
def plot_spike_EMG_xds(dataset, spikes, EMG, EMG_chs, offset, raw_flag, bin_size = 0, force = [], force_bin_size = []):
    if raw_flag == 1:
        try:
            bin_size = 1/dataset.EMG_fs
        except Exception:
            bin_size = stats.mode(np.diff(dataset.raw_EMG_time_frame))[0][0]
    else:
        bin_size = bin_size 
    if 'EMG' in dataset.EMG_names[0]:
        p_names = [each[4:] for each in dataset.EMG_names]
    else:
        p_names = dataset.EMG_names
    N = len(EMG_chs)
    spike_grid = 5
    grid = plt.GridSpec(N+spike_grid,1,wspace=0.5,hspace=0.2)
    main_ax = plt.subplot(grid[0:spike_grid,0])
    for i, spiketrain in enumerate(spikes):
        main_ax.plot(spiketrain - offset, np.ones_like(spiketrain) * i, ls='', marker='|', color = 'k', ms = 1)
    if force != []:
       force_ax = main_ax.twinx()
       t_force = np.arange(force.shape[0])*force_bin_size - offset
       force_ax.plot(t_force, force[:, 0], 'blue')
       force_ax.plot(t_force, force[:, 1], 'royalblue')
       force_ax.axis('off') 

    main_ax.axis('off')
    plt.xticks(color = 'w')
    plt.yticks([])
    #ylim_num = 1*np.max(p_emg[plot_start:plot_start+plot_len, :])
    x = np.arange(EMG.shape[0])*bin_size - offset
    ymin, ymax = np.min([EMG[:, EMG_chs[i]] for i in range(N)]), np.max([EMG[:, EMG_chs[i]] for i in range(N)])
    
    for i in range(N):
        ax0 = plt.subplot(grid[i+spike_grid,0], sharex = main_ax)
        p1 = EMG[:, EMG_chs[i]]
        #plt.yticks([])
        frame = plt.gca()
        frame.axes.get_yaxis().set_visible(False)
        ax0.spines['top'].set_visible(False)
        ax0.spines['right'].set_visible(False)
        ax0.spines['left'].set_visible(False)
        if i<N-1:
            plt.plot(x, p1, 'k')
            ax0.spines['bottom'].set_visible(False)
            plt.setp(ax0.get_xticklabels(),visible=False)
            ax0.tick_params(axis=u'both', which=u'both',length=0)
        if i == N-1:
            ax0.tick_params(axis=u'both', which=u'both',length=4)
            plt.setp(ax0.get_xticklabels(),visible=True)
            plt.plot(x, p1, 'k')
            ax0.set_xlabel('Time (s)', fontsize = 20)
            plt.tick_params(labelsize = 16)
            labels = ax0.get_xticklabels() + ax0.get_yticklabels()
            [label.set_fontname('Arial') for label in labels]
        plt.ylim(ymin, ymax)
        plt.text(x[-1], np.max(p1),'%s' %(p_names[EMG_chs[i]]),fontsize = 16, color = 'k',
                  verticalalignment="top",horizontalalignment="left")


def plot_EMGs(arr1, bin_width):
    T, N = arr1.shape[0], arr1.shape[1]
    t = np.arange(T)*bin_width
    for i in range(N):
        plt.subplot(N, 1, i+1)
        plt.plot(t, arr1[:, i], color = 'k', linewidth = 2)
        plt.ylim([np.min(arr1[:, i]), 1.2*np.max(arr1[:, i])])
    
def plot_EMGs_distribution(arr1, title_str):
    N = arr1.shape[1]  
    for i in range(N):
        plt.subplot(2, 4, i+1)
        ax = sns.histplot(arr1[:, i], kde=True, color = 'dimgray')
        ax.set_title(title_str[i])
        plt.tight_layout()
        # plt.xlim([-0.2, 2])
        
    
    
    
    
    
    
    
