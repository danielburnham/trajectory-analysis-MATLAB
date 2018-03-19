%% 4. this program will read in a set of values to make the ref subtracted traces start at zero

clear all; close all; clc;

%% paths
path = '/Volumes/DanielBurnham';

%% variables to change
filter_N = 348;                                                     % filter window size
date = 'yyyy/yyyy-mm-dd';                                           % date of experiment
filt = ['z_ref_sub_filt_filter_N_' num2str(filter_N) '.dat'];       % name of filtered data file
zl = 'zl_offsets.dat';                                              % name of file while zero offsets are stored
export_file_name = ['z_zl_filter_N_' num2str(filter_N) '.dat'];     % name of file in which to save data

%% read filtered data
z_ref_sub_filt = csvread([path '/' date '_analysis' '/' filt]);

%% read zero levels
zero_levels = csvread([path '/' date '_analysis' '/zl/' zl]);

%% use zero levels to get the traces to start at zero
z_zl = bsxfun(@minus,z_ref_sub_filt,zero_levels);

%% output the zero leveled z data ref subtracted and filtered for all beads
filenametosave = [path '/' date '_analysis' '/' export_file_name];
dlmwrite(filenametosave,z_zl,'newline','pc','precision','%.6f');