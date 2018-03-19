% this program takes the z reference subtracted data and filters it,
% outputs this new data and the zero level for all beads; for the constant
% zero level trace

clear all; close all; clc;

%% paths
path = '/Volumes/DanielBurnham';

%% variables to change
date = 'yyyy/yyyy-mm-dd';                   %date of experiment
z_ref_sub_data = 'zl_z_ref_sub.dat';        % name of file in which the z ref subtracted data is stored
export_file_name = 'zl_z_ref_sub_filt.dat'; % name of file in which to save data
zero_level_offsets = 'zl_offsets.dat';      % name of file in which to save zero level values
fps = 58;                                   % frame rate of raw data
filter_N = 348;                             % filter window size
dt = 1/fps;

%% read ref subtracted z data
z_ref_sub = csvread([path '/' date '_analysis' '/zl/' z_ref_sub_data]);

%% find number of beads
no_beads = size(z_ref_sub,2);

%% moving mean filter
for i = 1:no_beads  
    z_ref_sub_filt(:,i) = movmean(z_ref_sub(:,i),filter_N);
end

%% output ref subtracted and filtered z for all beads
filenametosave = [path '/' date '_analysis' '/zl/' export_file_name];
dlmwrite(filenametosave,z_ref_sub_filt,'newline','pc','precision','%.6f');

%% key difference is finding mean of zero level
zero_level = mean(z_ref_sub_filt,1);

%% output zero level value for all beads
filenametosave2 = [path '/' date '_analysis' '/zl/' zero_level_offsets];
dlmwrite(filenametosave2,zero_level,'newline','pc','precision','%.6f');