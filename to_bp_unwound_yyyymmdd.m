%% 5. this program will convert zero leveled traces to nm then to bp unwound

clear all; close all; clc;

%% paths
path = '/Volumes/DanielBurnham';

%% variables to change
filter_N = 348;                                                     % filter window size
date = 'yyyy/yyyy-mm-dd';                                           % date of experiment
z_zl = ['z_zl_filter_N_' num2str(filter_N) '.dat'];                 % name of file for zeroed data;
export_file_name = ['z_bp_filter_N_' num2str(filter_N) '.dat'];     % name of file to export filtered data to
LUT_step = 100;                                                     % LUT step size in nm
interface_factor = 0.88;                                             % actual focal shift due to mismatched refractive indices

%% read in z zero level data
z = csvread([path '/' date '_analysis' '/' z_zl]);

%% convert to nm, then bp
z_nm = z * LUT_step * interface_factor;                             % convert to nm
z_bp = (z_nm/((1.5088-0.9716)*1000))*2711;                          % convert to bp unwound (z_meas - z_ds / z_ss - z_ds)*2700 = (z_measured/((1.5088-0.9716)*1000)) * 2700

%% output the bp unwound data zero leveled data ref subtracted and filtered z for all beads
filenametosave = [path '/' date '_analysis' '/' export_file_name];
dlmwrite(filenametosave,z_bp,'newline','pc','precision','%.6f');