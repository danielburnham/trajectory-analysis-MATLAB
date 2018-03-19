%% this program takes the raw tracked data and outputs z traces only; for the constant zero level trace

clear all; close all; clc;

%% pre-assign for growing in loop
xyz = [];

%% paths
    read_path = '/Volumes/DanielBurnham';
    write_path = '/Volumes/DanielBurnham';

%% variables to change
date = 'yyyy/yyyy-mm-dd'; % date of experiment
tracked_data = 'raw.dat'; % name of tracked data file
export_file_name = 'zl_z_only_all_beads.dat'; % name of file in which to save data

%% read raw data
for j = 9:9 %folders that contain the force constant pre-unwinding data
    filename_zl = [read_path '/' date '/goagain' num2str(j) '/' tracked_data];
    xyz_raw = csvread(filename_zl);
    xyz = [xyz; xyz_raw];
end

%% find number of beads
no_beads = size(xyz,2)/3;

%% extract only z columns
z = xyz(:,3*linspace(1,no_beads,no_beads));

%% make a directory
[status, msg, msgID] = mkdir([write_path '/' date '_analysis' '/zl/']);

%% output only z for all beads
filenametosave = [write_path '/' date '_analysis' '/zl/' export_file_name];
dlmwrite(filenametosave,z,'newline','pc','precision','%.6f');