%% 6. this program will take the bp unwound traces, remove the non-unwinding parts and export as individual traces

clear all; close all; clc;

%% paths
path = '/Volumes/DanielBurnham';

%% variables to change
filter_N = 348;                                         % filter window size
fps = 58;
dt = 1/fps;
date = 'yyyy/yyyy-mm-dd';                               % date of experiment
z_bp = ['z_bp_filter_N_' num2str(filter_N) '.dat'];     % name of file for bp unwound data;
interval = 1;                                           % later if only every n^th frame tracked
st_fold_no = 9;
max_folder_no = 71;

%% read in z bp unwound data level data
bp_unwound = csvread([path '/' date '_analysis' '/' z_bp]);

%% beads that are for analysis
bead = [1 2 3 8 10 14 15 16 17 18 20 21 25 26 28 29 31 35 37 38 41 43 44 45 46];

%% folder where EVERY bead disappears including the non-analysed ones (i.e. we know it was in the previous folder at frame 3999)
% bead_dis = [37,8,7,84,60,4,25,84,84,4,8,11,14,7,27,23,84,84,84,5,5,36,36,84,84,84,84,84,84,84,84,41,33,44,39,40,4,14,7,84,84,8,24,84,84,84,36,33,84,71,10,84,84,84,22,41,3,32,84,84,84,6,25,13,12,9,13,20,32,5,5,4,14,17,84,19,29,61,84,84];
bead_dis = 71*ones(1,size(bp_unwound,2));

%% array of numbers that say where only the beads that are worthy of being tracked have disappeared
%  i.e. bead array and bdt now correlate with each other
bdt = bead_dis(bead);

%% the same as above except now it is the folder in which the bead is still there.
bst = bdt - 1;

%% check it makes sense
if size(bead) == size(bdt)
    % do nothing
else
    error('oops')
end

%% sets up start_points and end_points to extract the correct part of the trace for each trajectory
start_points = ones(1,length(bead));
end_points = floor(max_folder_no*4000*ones(1,size(bead,2))/interval);

for j = 1:size(bp_unwound,2)
    list_no = find(bead == j);
    end_points(list_no) = floor((bst(list_no)-st_fold_no+1)*4000/interval)l;
end

%% crops each trajectory as per the above start and end arrrays and outputs the individual trajectory files
for i = 1:size(bead,2)
    disp(i);
    if end_points(i) < 1
        bp_unwound_cut = bp_unwound(start_points(i):start_points(i)+1,bead(i));
    else
        bp_unwound_cut = bp_unwound(start_points(i):end_points(i),bead(i));
    end
    
    %% time vector for each molecule
    time = linspace(0,size(bp_unwound_cut,1)-1,size(bp_unwound_cut,1)).'*dt*interval;

    %% clears out NaNs for fitting of line
    validdata1 = ~isnan(time);
    validdata2 = ~isnan(bp_unwound_cut);
    validdataBoth = validdata1 & validdata2;
    time_nanless = time(validdataBoth) ;
    bp_co_nanless = bp_unwound_cut(validdataBoth);
    
    figure(2)
    plot(time_nanless,bp_co_nanless)
    hold on
    
    %% export nanless time
    filenametosave1 = [path '/' date '_analysis' '/' 'time' '_filter_N_' num2str(filter_N) '_mol_' num2str(bead(i)) '.dat'];
    dlmwrite(filenametosave1,time_nanless,'newline','pc','precision','%.6f');
    
    %% export nanless trace
    filenametosave2 = [path '/' date '_analysis' '/' 'bp_final' '_filter_N_' num2str(filter_N) '_mol_' num2str(bead(i)) '.dat'];
    dlmwrite(filenametosave2,bp_co_nanless,'newline','pc','precision','%.6f');
end