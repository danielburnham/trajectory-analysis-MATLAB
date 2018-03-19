%% 8. this program will calcualte all passage times and plot them

clear all; close all; clc;

%% paths
path = '/Volumes/DanielBurnham';

%% variables to change
date = 'yyyy/yyyy-mm-dd';               % date of experiment
time_file = ['time' '_filter_N_'];
z_file = ['bp_final' '_filter_N_'];     % name of file for individual traces of bp unwound data;
fps = 58;
dt = 1/fps;
filter_N = 348;
interval = 1;

%% choose first passage interval (bp)
int = 20;

%% set up empty array
all_PASSTIME = [];

%% read in beads to analyse
beads_ta = csvread([path '/' date '_analysis' '/' 'beads_to_analyse' '_filter_N_' num2str(filter_N) '.dat']);

%% t for sub plot indices
t=0;

%% figure for all data
figure(1)

%% for loop over all beads that are in "beads to analyse"
for j = beads_ta
    disp(j)
    
    PASSTIME = [];
    tau = [];
    
    %% read in time data
    time = csvread([path '/' date '_analysis' '/' time_file num2str(filter_N) '_mol_' num2str(j) '.dat']);

    %% read in z_bp data
    z_bp = csvread([path  '/' date '_analysis' '/' z_file num2str(filter_N) '_mol_' num2str(j) '.dat']);

    %% plot
    t=t+1;
    subplot(9,9,t)
    plot(time,z_bp);
    hold on
    
    start_height = z_bp(1);                         % starting passage height
    
    gapdan = start_height + int;                    % next passage height
    
    %% calculate linear m and c from y=mx+c for each pair of points
    for q = 1:length(z_bp)-1
        m_interp(q) = (z_bp(q+1)-z_bp(q))/((time(q+1))-(time(q)));
        c_interp(q) = z_bp(q) - (m_interp(q)*(time(q)));
    end
    
    %% plot horizontal passage intervals and vertical passage crossings
    hoz = start_height;
    k = 1;
    plot([0,6000],[hoz,hoz],'g-')
    
    for i = 1:length(z_bp)-1        
        while hoz <= z_bp(i+1)           
            k = k + 1;
            tau(k-1) = (hoz-c_interp(i))/m_interp(i);
            
            plot([0,6000],[hoz,hoz],'g-')
            plot([tau(k-1),tau(k-1)],[0,6000],'k-')
            
            hoz = hoz + int;
        end      
    end
    
    axis([0 max(time) 0 max(z_bp)])
    
    for g = 1:length(tau)-1
        PASSTIME(g) = tau(g+1) - tau(g);
    end
    
    all_PASSTIME = [all_PASSTIME; PASSTIME.'];
    
    
end

%% output all passage times
filenametosave = [path '/' date '_analysis' '/' 'all_PASSTIME' '_filterN_' num2str(filter_N) '.dat'];
dlmwrite(filenametosave,all_PASSTIME,'newline','pc','precision','%.6f');

%% plot hisotgram of passage times
% make bins that are equally spaced in log space
logbins = logspace(-1,4,28);

% count things in each of these log spaced bins
[logcounts,logedges] = histcounts(all_PASSTIME,logbins);

% what is the width of these log spacd bins?
for i = 2:length(logbins)
    logwidth(i-1) = logbins(i) - logbins(i-1);
end

% the bars are normalised by total counts and width to make a probability density function
logbars = ((logcounts./sum(logcounts))./logwidth);

figure(2)
loglog(logedges(1:end-1),logbars,'bo','MarkerSize',10,'LineWidth',1.5)
axis([0.0001 1E4 1E-8 10])