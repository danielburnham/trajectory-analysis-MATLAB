%% 2. this program takes the z only data and subtracts the reference from it all

clear all; close all; clc;

%% paths
path = '/Volumes/DanielBurnham';

%% variables to change
date = 'yyyy/yyyy-mm-dd';                       % date of experiment
z_only_tracked_data = 'z_only_all_beads.dat';   % name of file in which the z only data is stored
export_file_name = 'z_ref_sub.dat';             % name of file in which to save data
ref = 44;                                       % ref bead

%% read z only tracked data
z_raw = csvread([path '/' date '_analysis' '/' z_only_tracked_data]);

%% subtract ref from all beads
z_ref_sub = bsxfun(@minus,z_raw,z_raw(:,ref));

%% output ref subtract z for all beads
filenametosave = [path '/' date '_analysis' '/' export_file_name];
dlmwrite(filenametosave,z_ref_sub,'newline','pc','precision','%.6f');