%% 1. this program takes the raw xyz tracked data and outputs z data only

clear all; close all; clc;

%% pre-assign fro growing in loop
xyz = [];

%% paths
read_path = '/Volumes/DanielBurnham';
write_path = '/Volumes/DanielBurnham';

%% variables to change
date = 'yyyy/yyyy-mm-dd';                   % date of experiment
tracked_data = 'raw.dat';                   % name of tracked data file
export_file_name = 'z_only_all_beads.dat';  % name of file in which to save data

%% read tracked data from all folders
for j = 9:70
    disp(j);
    filename = [read_path '/' date '/goagain' num2str(j) '/' tracked_data];
    xyz_raw = csvread(filename);
    xyz = [xyz; xyz_raw];
end

%% find number of beads
no_beads = size(xyz,2)/3;

%% extract only z columns
z = xyz(:,3*linspace(1,no_beads,no_beads));

%% output only z for all beads
filenametosave = [write_path '/' date '_analysis' '/' export_file_name];
dlmwrite(filenametosave,z,'newline','pc','precision','%.6f');