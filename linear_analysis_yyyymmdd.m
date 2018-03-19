%% 7. performs a linear fit to the data

clear all; close all; clc;

%% paths
addpath(genpath('plotSpread'))
path = '/Volumes/DanielBurnham';

%% variables to change
filter_N = 348;                                         % filter window size
fps = 58;
dt = 1/fps;
date = 'yyyy/yyyy-mm-dd';                               % date of experiment
partial = ['_filter_N_' num2str(filter_N) '_mol_'];     % name of file for bp unwound data;
interval = 1;
st_fold_no = 9;
ref = 44;
gradient_thresh = 0.001;


%% beads that are for analysis
bead = [1 2 3 8 10 14 15 16 17 18 20 21 25 26 28 29 31 35 37 38 41 43 44 45 46];

t = 0
for i = bead
    t = t+1
    
    %% read in bp unwound individual trajectories
    time = csvread([path '/' date '_analysis' '/' 'time' partial num2str(i) '.dat']);
    bp_unwound = csvread([path '/' date '_analysis' '/' 'bp_final' partial num2str(i) '.dat']);
    
    %% use mldivide ('\') to perform linear fit through origin
    grad = time(:)\bp_unwound(:)
    
    %% evalulate linear fit for plotting
    x_for_fit = [time(1) time(end)];
    y_eval = grad*x_for_fit
    
    %% for plottting many trajectories
    if i <= 42
        f1 = figure(1);
        subplot(7,6,i)
    elseif i >42 && i <= 84
        f11 = figure(11);
        subplot(7,6,i-42)
    else
        f111 = figure(111);
        subplot(7,6,i-84)
    end
    
    plot(time,bp_unwound)
    drawnow
    hold on
    plot(x_for_fit,y_eval,'r-','LineWidth',1.2)
    
    %% make 'm's zero if total time less than one second so they get ignored
    if max(time) < 1
        m(t) = 0
    else
        m(t) = grad;
    end
    
    title(['m=' num2str(m(t),'%0.2f') ',' 'b# =' num2str(i)]);
    
end

% sort out which beads are unwinding
% remove ref beads
[ynref,locref] = ismember(ref,bead);
bead(locref) = [];
m(locref) = [];

% sort out which gradients are unwinding
beads_to_analyse = bead(m > gradient_thresh);
[ynm,locm] = ismember(beads_to_analyse,bead);
m_above_thresh = m(locm);

if length(m_above_thresh) < 1
    %do nothing
else
    %% find mean and SEM
    mean_m = mean(m_above_thresh)
    std_err_mean = std(m_above_thresh)/sqrt(length(m_above_thresh)-1)
    
    %% plot data
    f2 = figure(2);
    p1 = plotSpread(m_above_thresh.','distributionMarkers','o','binWidth',0.1);
    
    oldpos = get(f2,'Position');
    newpos = [oldpos(1),oldpos(2)-200,250,560];
    set(f2,'Position',newpos)
    
    t = allchild(gca);
    set(t,'MarkerSize',10,'LineWidth',1.2);
    
    box on
    hold on
    b1 = boxplot(m_above_thresh,'Widths',1.0,'Color','k','Whisker',10)
    
    for ib = 1:6
        set(b1(ib,:),'LineWidth',1.2)
    end
    
    set(gca,'XTickLabel','')
    ylabel('Velocity (bps^{-1})','fontsize',14)
    set(gca,'fontsize',12,'linewidth',1.2)
    
    uistack(b1,'bottom')
    to_add = strcat([num2str(mean_m,'%0.2f') ' \pm ' num2str(std_err_mean,'%0.1g') ' bps^{-1}'])
    text(0.15,0.525,to_add,'fontsize',12)
end

%% output list of beads_to_analyse in fptd
filenametosave = [path '/' date '_analysis' '/' 'beads_to_analyse' '_filter_N_' num2str(filter_N) '.dat'];
dlmwrite(filenametosave,beads_to_analyse,'newline','pc','precision','%u');

%% output list of gradients
filenametosave2 = [path '/' date '_analysis' '/' 'gradients' '_filter_N_' num2str(filter_N) '.dat'];
dlmwrite(filenametosave2,m_above_thresh,'newline','pc','precision','%.6f');