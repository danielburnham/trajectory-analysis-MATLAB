%% 3. this program takes the z reference subtracted data and filters it

clear all; close all; clc;

%% paths
path = '/Volumes/DanielBurnham';

%% variables to change
fps = 58;                                                                   % frame rate of raw data
filter_N = 348;                                                             % filter window size
dt = 1/fps;
date = 'yyyy/yyyy-mm-dd';                                                   % date of experiment
z_ref_sub_data = 'z_ref_sub.dat';                                           % name of file in which the z ref subtracted data is stored
export_file_name = ['z_ref_sub_filt_filter_N_' num2str(filter_N) '.dat'];   % name of file in which to save data

%% read ref subtracted z data
z_ref_sub = csvread([path '/' date '_analysis' '/' z_ref_sub_data]);

%% find number of beads
no_beads = size(z_ref_sub,2);

for i = 1:no_beads
    z_ref_sub_filt(:,i) = movmean(z_ref_sub(:,i),filter_N);
end

%% output ref subtracted and filtered z for all beads
filenametosave = [path '/' date '_analysis' '/' export_file_name];
dlmwrite(filenametosave,z_ref_sub_filt,'newline','pc','precision','%.6f');