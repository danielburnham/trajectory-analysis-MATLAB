function [A,B,lb_for_fit_truncated,combined_fit_1_double,combined_fit_2_double,combined_fit_3_double,logedges,logbars,bic] = mle_fptd_yyyymmdd(laces_double,laces_triple,se_fr,de_fr,e3_fr,e4_fr,no_boots)
% close all; clear all; clc;
% laces_double: bootstrap for 2 exponentials? 0 = false, 1 = true
% laces_triple: bootstrap for 3 exponentials? 0 = false, 1 = true
% se_fr: fit with 1 exponentials
% de_fr: fit with 2 exponentials
% e3_fr: fit with 3 exponentials
% e4_fr: fit with 4 exponentials
% no_boots: number of bootstraps
%% this function fits different models to the data and bootstraps

colors1 = lbmap(50,'RedBlue');
colors2 = lbmap(50,'BrownBlue');
colors3 = lbmap(50,'BlueGray');

%% bootstrap yes or no
% 0 = no
% 1 = yes
% laces = 1;

%% single exponential fitting rouitne
% 0 = off
% 1 = on
% se_fr = 0;

%% double exponential fitting routine
% 0 = off
% 1 = on
% de_fr = 1;

%% 3 exponential fitting routine
% 0 = off
% 1 = on
% e3_fr = 0;

%% settings
filter_N = 348;
interval_bp = 20;
gen1 = 10000;
gen2 = 10000;
gen3 = 10000;
tol = 1E-12;
% no_boots = 3;

%% paths
% path = '/Volumes/DanielBurnham';
read_path = '/Volumes/DanielBurnham';
date = 'yyyy/yyyy-mm-dd';
filename = [read_path '/' date '_analysis/' 'all_PASSTIME' '_filterN_' num2str(filter_N) '.dat']

%% get passage times
all_PASSTIME = csvread(filename);
pass_times = all_PASSTIME.';
lb_for_fit = logspace(-4,4,1000);
lb_for_fit_truncated = logspace(log10(min(pass_times))-0.5,log10(max(pass_times)),1000);

%% make probability distribution
% make bins that are equally spaced in log space
logbins = logspace(-2,4,32);
% count things in each of these log spaced bins
[logcounts,logedges] = histcounts(all_PASSTIME,logbins);

% what is the width of these log spacd bins?
for i = 2:length(logbins)
    logwidth(i-1) = logbins(i) - logbins(i-1);
end

% the bars are normalised by total counts and width to make a probability density function
logbars = ((logcounts./sum(logcounts))./logwidth);
figure(211)
loglog(logedges(1:end-1),logbars,'bo','MarkerSize',10,'LineWidth',1.5)
hold on
axis([0.01 1E4 1E-8 10])
drawnow

%% single exponential fit for combined function passage time
oop1 = optimoptions('fmincon','StepTolerance',1E-21,'ConstraintTolerance',1E-21,'OptimalityTolerance',1E-21);

% set up linear equality constraints
acon = [0 0 0 1 1];
bcon = 1;

if se_fr == 0
    % do nothing
    bic_one = 0;
    A = zeros(1,6);
elseif se_fr == 1
    % set up low and high constraints
    rf_low = 10;        %A(1)
    rb_low = 10;        %A(2)
    k1_low = 0.0001;     %A(3)
    p1_low = 0;         %A(4)
    p2_low = 0;         %A(5)
    
    rf_high = 600;
    rb_high = 600;
    k1_high = 0.1;
    p1_high = 1;
    p2_high = 1;
    
    % use genetic algorithm to perform MLE
    [A fval exitflag] = fmincon(@(A) -LLexp_and_LLRW(A,length(pass_times),pass_times,interval_bp),[110 80 0.01 0.5 0.5],[],[],acon,bcon,[rf_low rb_low k1_low p1_low p2_low],[rf_high rb_high k1_high p1_high p2_high],[],oop1)
    
    % evaluate model using estiamted parameters and lb_for_fit times
    first_fit = (interval_bp./lb_for_fit_truncated)*((A(1)/A(2))^(interval_bp/2)).*exp(((-((A(1) + A(2))/1))+(( (A(1) + A(2))*sqrt(1-(((A(1) - A(2))/(A(1) + A(2)))^2)) ))).*lb_for_fit_truncated);
    second_fit = besseli(interval_bp,  ((A(1) + A(2))/1)*sqrt(1-(((A(1) - A(2))/(A(1) + A(2)))^2))*lb_for_fit_truncated,1    );
    third_fit = A(3)*exp(-A(3).*lb_for_fit_truncated);
    
    combined_fit_1 = A(4).*(first_fit.*second_fit);     % RW
    combined_fit_2 = A(5).*third_fit;                   % pauses
    
    % calculate Bayesian information criterion
    bic_one = LLexp_and_LLRW(A,length(pass_times),pass_times,interval_bp);
    
end

%% double exponential fit for combined function passage time
oop2 = optimoptions('fmincon','StepTolerance',1E-21,'ConstraintTolerance',1E-21,'OptimalityTolerance',1E-21);

% set up linear equality constraints
acon = [0 0 0 0 1 1 1];
bcon = 1;

% set up linear inequality constraints
aincon = [0 0 -1 1 0 0 0; 0 0 0 0 0 0 0];
bincon = [0; 0];
% evaluates the matrix product A*x as if x is transposed (A*x').

if de_fr == 0
    % do nothing
    bic_two = 0;
    B = zeros(1,9);
elseif de_fr == 1;
    % set up low and high constraints
    rf_low_dub = 10;      % B(1)
    rb_low_dub = 10;      % B(2)
    k1_low_dub = 0.0001;	% B(3)
    k2_low_dub = 0.0001;	% B(4)
    p1_low_dub = 0;         % B(5)
    p2_low_dub = 0;         % B(6)
    p3_low_dub = 0;         % B(7)
    
    rf_high_dub = 600;
    rb_high_dub = 600;
    k1_high_dub = 0.1;
    k2_high_dub = 0.1;
    p1_high_dub = 1;
    p2_high_dub = 1;
    p3_high_dub = 1;
    
    % use fmincon to perform MLE
    [B] = fmincon(@(B) -LLdoubleexp_and_LLRW(B,length(pass_times),pass_times,interval_bp),[110 80 0.01 0.001 0.5 0.5 0.5],aincon,bincon,acon,bcon,[rf_low_dub rb_low_dub k1_low_dub k2_low_dub p1_low_dub p2_low_dub p3_low_dub],[rf_high_dub rb_high_dub k1_high_dub k2_high_dub p1_high_dub p2_high_dub p3_high_dub],[],oop2);
    
    first_fit_double = (interval_bp./lb_for_fit_truncated)*((B(1)/B(2))^(interval_bp/2)).*exp(((-((B(1) + B(2))/1))+(( (B(1) + B(2))*sqrt(1-(((B(1) - B(2))/(B(1) + B(2)))^2)) ))).*lb_for_fit_truncated);
    second_fit_double = besseli(interval_bp, ((B(1) + B(2))/1)*sqrt(1-(((B(1) - B(2))/(B(1) + B(2)))^2))*lb_for_fit_truncated,1    );
    third_fit_double = B(3)*exp(-B(3).*lb_for_fit_truncated);
    fourth_fit_double = B(4)*exp(-B(4).*lb_for_fit_truncated);
    
    combined_fit_1_double = B(5).*(first_fit_double.*second_fit_double);
    combined_fit_2_double = (B(6)).*third_fit_double;
    combined_fit_3_double = B(7).*fourth_fit_double;
    
    bic_two = LLdoubleexp_and_LLRW(B,length(pass_times),pass_times,interval_bp);
    actual_params = [read_path '/' date '_analysis/' '/' 'actual_params_2.csv']
    dlmwrite(actual_params,B.','newline','pc','precision','%.8f')
    
end

%% 3 exponential fit for combined function passage time
oop3 = optimoptions('fmincon','StepTolerance',1E-21,'ConstraintTolerance',1E-21,'OptimalityTolerance',1E-21);

% set up linear equality constraints
acon = [0 0 0 0 0 1 1 1 1];
bcon = 1;

% set up linear inequality constraints
aincon = [0 0 0 -1 1 0 0 0 0; 0 0 -1 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0];
bincon = [0; 0; 0; 0];

if e3_fr == 0
    % do nothing
    bic_three = 0;
    C = zeros(1,12);
elseif e3_fr == 1
    % set up low and high constraints
    rf_low = 10;
    rb_low = 10;
    k1_low = 0.0001;
    k2_low = 0.0001;
    k3_low = 0.0001;
    p1_low = 0;
    p2_low = 0;
    p3_low = 0;
    p4_low = 0;
    
    rf_high = 600;
    rb_high = 600;
    k1_high = 0.1;
    k2_high = 0.1;
    k3_high = 0.1;
    p1_high = 1;
    p2_high = 1;
    p3_high = 1;
    p4_high = 1;
    
    % use fmincon to perform MLE
    [C] = fmincon(@(B3) -LL_3_exp_and_LLRW(B3,length(pass_times),pass_times,interval_bp),[110 80 0.01 0.001 0.0001 0.5 0.5 0.5 0.5],aincon,bincon,acon,bcon,[rf_low rb_low k1_low k2_low k3_low p1_low p2_low p3_low p4_low],[rf_high rb_high k1_high k2_high k3_high p1_high p2_high p3_high p4_high],[],oop3);
    
    first_fit_3 = (interval_bp./lb_for_fit_truncated)*((C(1)/C(2))^(interval_bp/2)).*exp(((-((C(1) + C(2))/1))+(( (C(1) + C(2))*sqrt(1-(((C(1) - C(2))/(C(1) + C(2)))^2)) ))).*lb_for_fit_truncated);
    second_fit_3 = besseli(interval_bp, ((C(1) + C(2))/1)*sqrt(1-(((C(1) - C(2))/(C(1) + C(2)))^2))*lb_for_fit_truncated,1    );
    
    third_fit_3 = C(3)*exp(-C(3).*lb_for_fit_truncated);
    fourth_fit_3 = C(4)*exp(-C(4).*lb_for_fit_truncated);
    fifth_fit_3 = C(5)*exp(-C(5).*lb_for_fit_truncated);
    
    combined_fit_1_3 = C(6).*(first_fit_3.*second_fit_3);
    combined_fit_2_3 = C(7).*third_fit_3;
    combined_fit_3_3 = C(8).*fourth_fit_3;
    combined_fit_4_3 = C(9).*fifth_fit_3;
    
    bic_three = LL_3_exp_and_LLRW(C,length(pass_times),pass_times,interval_bp);
    
end

%% 4 exponential fit for combined function passage time
oop4 = optimoptions('fmincon','StepTolerance',1E-21,'ConstraintTolerance',1E-21,'OptimalityTolerance',1E-21);

% set up linear equality constraints
acon = [0 0 0 0 0 0 1 1 1 1 1];
bcon = 1;

% set up linear inequality constraints
aincon = [0 0 -1 1 0 0 0 0 0 0 0; 0 0 0 -1 1 0 0 0 0 0 0; 0 0 0 0 -1 1 0 0 0 0 0];
bincon = [0; 0; 0];

if e4_fr == 0
    % do nothing
    bic_four = 0;
    D = zeros(1,12);
elseif e4_fr == 1
    % set up low and high constraints
    rf_low = 10;
    rb_low = 10;
    k1_low = 0.0001;
    k2_low = 0.0001;
    k3_low = 0.0001;
    k4_low = 0.0001;
    p1_low = 0;
    p2_low = 0;
    p3_low = 0;
    p4_low = 0;
    p5_low = 0;
    
    rf_high = 600;
    rb_high = 600;
    k1_high = 0.1;
    k2_high = 0.1;
    k3_high = 0.1;
    k4_high = 0.1;
    p1_high = 1;
    p2_high = 1;
    p3_high = 1;
    p4_high = 1;
    p5_high = 1;
    
    % use fmincon to perform MLE
    [D] = fmincon(@(B3) -LL_4_exp_and_LLRW(B3,length(pass_times),pass_times,interval_bp),[110 80 0.01 0.001 0.0001 0.00001 0.5 0.5 0.5 0.5 0.5],aincon,bincon,acon,bcon,[rf_low rb_low k1_low k2_low k3_low k4_low p1_low p2_low p3_low p4_low p5_low],[rf_high rb_high k1_high k2_high k3_high k4_high p1_high p2_high p3_high p4_high p5_high],[],oop4);
    
    first_fit_4 = (interval_bp./lb_for_fit_truncated)*((D(1)/D(2))^(interval_bp/2)).*exp(((-((D(1) + D(2))/1))+(( (D(1) + D(2))*sqrt(1-(((D(1) - D(2))/(D(1) + D(2)))^2)) ))).*lb_for_fit_truncated);
    second_fit_4 = besseli(interval_bp, ((D(1) + D(2))/1)*sqrt(1-(((D(1) - D(2))/(D(1) + D(2)))^2))*lb_for_fit_truncated,1    );
    
    third_fit_4 = D(3)*exp(-D(3).*lb_for_fit_truncated);
    fourth_fit_4 = D(4)*exp(-D(4).*lb_for_fit_truncated);
    fifth_fit_4 = D(5)*exp(-D(5).*lb_for_fit_truncated);
    sixth_fit_4 = D(6)*exp(-D(6).*lb_for_fit_truncated);
    
    combined_fit_1_4 = D(7).*(first_fit_4.*second_fit_4);
    combined_fit_2_4 = D(8).*third_fit_4;
    combined_fit_3_4 = D(9).*fourth_fit_4;
    combined_fit_4_4 = D(10).*fifth_fit_4;
    combined_fit_5_4 = D(11).*sixth_fit_4;
    
    bic_four = LL_4_exp_and_LLRW(D,length(pass_times),pass_times,interval_bp);
    
end

%% bootstrapping
mean_point = zeros(size(logbars,1),1);
std_point = zeros(size(logbars,1),1);
mu_param = zeros(size(C,2),1);
std_param = zeros(size(C,2),1);

%% bootstrap passage times and fit for double exponential
if laces_double == 0
    % do nothing
    %     mean_point = zeros(size(logbars,1),1);
    %     std_point = zeros(size(logbars,1),1);
    %     mu_param = zeros(size(B,2),1);
    %     std_param = zeros(size(B,2),1);
elseif laces_double == 1
    [~,bootsam] = bootstrp(no_boots,[],pass_times);
    figure(19)
    
    for d = 1:no_boots
        d
        % set up linear equality constraints
        acon = [0 0 0 0 1 1 1];
        bcon = 1;
        
        % set up linear inequality constraints
        aincon = [0 0 -1 1 0 0 0; 0 0 0 0 0 0 0];
        bincon = [0 0];
        
        % perform bootstrap
        bootstrap_no = d
        matlab_boot(:,d) = pass_times(bootsam(:,d));
        
        % use fmincon to perform MLE
        [BBB] = fmincon(@(BBB) -LLdoubleexp_and_LLRW(BBB,length(matlab_boot(:,d)),matlab_boot(:,d).',interval_bp),[110 80 0.01 0.001 0.5 0.5 0.5],aincon,bincon,acon,bcon,[rf_low_dub rb_low_dub k1_low_dub k2_low_dub p1_low_dub p2_low_dub p3_low_dub],[rf_high_dub rb_high_dub k1_high_dub k2_high_dub p1_high_dub p2_high_dub p3_high_dub],[],oop2);
        
        BB(:,d) = BBB;
        
        % for plotting all the bootstraps
        first_fit_double = (interval_bp./lb_for_fit_truncated)*((BB(1,d)/BB(2,d))^(interval_bp/2)).*exp(((-((BB(1,d) + BB(2,d))/1))+(( (BB(1,d) + BB(2,d))*sqrt(1-(((BB(1,d) - BB(2,d))/(BB(1,d) + BB(2,d)))^2)) ))).*lb_for_fit_truncated);
        second_fit_double = besseli(interval_bp, ((BB(1,d) + BB(2,d))/1)*sqrt(1-(((BB(1,d) - BB(2,d))/(BB(1,d) + BB(2,d)))^2))*lb_for_fit_truncated,1    );
        third_fit_double = BB(3)*exp(-BB(3,d).*lb_for_fit_truncated);
        fourth_fit_double = BB(4)*exp(-BB(4,d).*lb_for_fit_truncated);
        
        combined_fit_1_doubleBB(:,d) = BB(5,d).*(first_fit_double.*second_fit_double);
        combined_fit_2_doubleBB(:,d) = (BB(6,d)).*third_fit_double;
        combined_fit_3_doubleBB(:,d) = BB(7,d).*fourth_fit_double;
        
        [logcounts,logedges] = histcounts(matlab_boot(:,d),logbins);
        logcounts_test(:,d) = logcounts.';
        
        for i = 2:length(logbins)
            logwidth_boot(i-1) = logbins(i) - logbins(i-1);
        end
        
        to_sum = logcounts_test(:,d);
        
        logbars_boot(:,d) = ((logcounts_test(:,d)./sum(to_sum))./logwidth_boot.');
        
    end
    % calculate the histogram error bars
    mean_point = mean(logbars_boot,2);
    std_point = std(logbars_boot,0,2);
    
    % calculate the errors on the parameters
    mu_param = mean(BB,2);
    std_param = std(BB,0,2);
    
    x = logedges(1:end-1).';
    y = logbars.';
    eb = std_point;
    param_err = std_param;
    
    % save true experimental histogram
    filename1 = [read_path '/' date '_analysis/' '/' 'xy.csv']
    dlmwrite(filename1,[x y],'newline','pc','precision','%.8f');
    
    % save all bootstrap bars
    filename2 = [read_path '/' date '_analysis/'  '/' 'logbars_boot.csv']
    dlmwrite(filename2,logbars_boot,'newline','pc','precision','%.8f');
    
    % save all bootstrap parameter values
    filename3 = [read_path '/' date '_analysis/'  '/' 'bootstrap_params.csv']
    dlmwrite(filename3,BB,'newline','pc','precision','%.8f');
    
    % save bootstrap bars standard deviation
    filename4 = [read_path '/' date '_analysis/'  '/' 'bootstrap_stdev.csv']
    dlmwrite(filename4,std_point,'newline','pc','precision','%.8f');
    
    % save bootstrap parameter standard deviation
    filename5 = [read_path '/' date '_analysis/'  '/' 'bootstrap_stdev.csv']
    dlmwrite(filename5,std_param,'newline','pc','precision','%.8f')
    
    % save bootstrapped passage times
    filename6 = [read_path '/' date '_analysis/'  '/' 'bootstrap_passage_times.csv']
    dlmwrite(filename6,matlab_boot,'newline','pc','precision','%.8f')
end

%% bootstrap passage times and fit for triple exponential
if laces_triple == 0
    % do nothing
    %     mean_point = zeros(size(logbars,1),1);
    %     std_point = zeros(size(logbars,1),1);
    %     mu_param = zeros(size(C,2),1);
    %     std_param = zeros(size(C,2),1);
elseif laces_triple == 1
    [~,bootsam] = bootstrp(no_boots,[],pass_times);
    figure(20)
    
    for d = 1:no_boots
        d
        % set up linear equality constraints
        acon = [0 0 0 0 0 1 1 1 1];
        bcon = 1;
        
        % set up linear inequality constraints
        aincon = [0 0 0 -1 1 0 0 0 0; 0 0 -1 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0];
        bincon = [0; 0; 0; 0];
        
        % perform bootstrap
        bootstrap_no = d
        matlab_boot(:,d) = pass_times(bootsam(:,d));
        
        % use fmincon to perform MLE
        [CCC] = fmincon(@(CCC) -LL_3_exp_and_LLRW(CCC,length(matlab_boot(:,d)),matlab_boot(:,d),interval_bp),[110 80 0.01 0.001 0.0001 0.5 0.5 0.5 0.5],aincon,bincon,acon,bcon,[rf_low rb_low k1_low k2_low k3_low p1_low p2_low p3_low p4_low],[rf_high rb_high k1_high k2_high k3_high p1_high p2_high p3_high p4_high],[],oop3);
        
        CC(:,d) = CCC;
        
        
        
        first_fit_3CC = (interval_bp./lb_for_fit_truncated)*((CC(1,d)/CC(2,d))^(interval_bp/2)).*exp(((-((CC(1,d) + CC(2,d))/1))+(( (CC(1,d) + CC(2,d))*sqrt(1-(((CC(1,d) - CC(2,d))/(CC(1,d) + CC(2,d)))^2)) ))).*lb_for_fit_truncated);
        second_fit_3 = besseli(interval_bp, ((CC(1,d) + CC(2,d))/1)*sqrt(1-(((CC(1,d) - CC(2,d))/(CC(1,d) + CC(2,d)))^2))*lb_for_fit_truncated,1    );
        
        third_fit_3CC = CC(3,d)*exp(-CC(3,d).*lb_for_fit_truncated);
        fourth_fit_3CC = CC(4,d)*exp(-CC(4,d).*lb_for_fit_truncated);
        fifth_fit_3CC = CC(5,d)*exp(-CC(5,d).*lb_for_fit_truncated);
        
        combined_fit_1_3CC = CC(6,d).*(first_fit_3.*second_fit_3);
        combined_fit_2_3CC = CC(7,d).*third_fit_3;
        combined_fit_3_3CC = CC(8,d).*fourth_fit_3;
        combined_fit_4_3CC = CC(9,d).*fifth_fit_3;
        
        [logcounts,logedges] = histcounts(matlab_boot(:,d),logbins);
        logcounts_test(:,d) = logcounts.';
        
        for i = 2:length(logbins)
            logwidth_boot(i-1) = logbins(i) - logbins(i-1);
        end
        
        to_sum = logcounts_test(:,d);
        
        logbars_boot(:,d) = ((logcounts_test(:,d)./sum(to_sum))./logwidth_boot.');
        
    end
    % calculate the histogram error bars
    mean_point = mean(logbars_boot,2);
    std_point = std(logbars_boot,0,2);
    
    % calculate the errors on the parameters
    mu_param = mean(CC,2);
    std_param = std(CC,0,2);
    
    x = logedges(1:end-1).';
    y = logbars.';
    eb = std_point;
    param_err = std_param;
    
    % save true experimental histogram
    filename1 = [read_path '/' date '_analysis/'  '/' 'xy3.csv']
    dlmwrite(filename1,[x y],'newline','pc','precision','%.8f');
    
    % save all bootstrap bars
    filename2 = [read_path '/' date '_analysis/'  '/' 'logbars_boot3.csv']
    dlmwrite(filename2,logbars_boot,'newline','pc','precision','%.8f');
    
    % save all bootstrap parameter values
    filename3 = [read_path '/' date '_analysis/'  '/' 'bootstrap_params3.csv']
    dlmwrite(filename3,CC,'newline','pc','precision','%.8f');
    
    % save bootstrap bars standard deviation
    filename4 = [read_path '/' date '_analysis/'  '/' 'bootstrap_stdev3.csv']
    dlmwrite(filename4,std_point,'newline','pc','precision','%.8f');
    
    % save bootstrap parameter standard deviation
    filename5 = [read_path '/' date '_analysis/'  '/' 'bootstrap_stdev3.csv']
    dlmwrite(filename5,std_param,'newline','pc','precision','%.8f')
    
    % save bootstrapped passage times
    filename6 = [read_path '/' date '_analysis/'  '/' 'bootstrap_passage_times3.csv']
    dlmwrite(filename6,matlab_boot,'newline','pc','precision','%.8f')
end


%% plot
f5 = figure(5);
axes1 = axes('Parent',f5);
hold(axes1,'on');
if laces_double == 0 & laces_triple == 0
    % plot without error bars
    h = loglog(logedges(1:end-1),logbars,'o','color',colors2(42,:),'MarkerFaceColor','w','MarkerSize',12,'LineWidth',1.4);
elseif laces_double == 1 | laces_triple == 1
    % plot with error bars
    h = errorbar(logedges(1:end-1),logbars,std_point,'o','color',colors2(42,:),'MarkerFaceColor','w','MarkerSize',12,'LineWidth',1.4);
    h.LData = h.YData - max(eps,h.YData - h.LData);
end

hold on

% plot single exponential version
if se_fr == 0;
    % do nothing
else
    loglog(lb_for_fit_truncated,combined_fit_1+combined_fit_2,'k-','Linewidth',1)
end

% plot double exponential version
if de_fr == 0;
    % do nothing
else
    loglog(lb_for_fit_truncated,combined_fit_1_double+combined_fit_2_double+combined_fit_3_double,'m-','Linewidth',1)
end

% plot triple exponential version
if e3_fr == 0;
    % do nothing
else
    loglog(lb_for_fit_truncated,combined_fit_1_3+combined_fit_2_3+combined_fit_3_3+combined_fit_4_3,'b-','Linewidth',1)
end

% plot triple exponential version
if e4_fr == 0;
    % do nothing
else
    loglog(lb_for_fit_truncated,combined_fit_1_4+combined_fit_2_4+combined_fit_3_4+combined_fit_4_4+combined_fit_5_4,'g-','Linewidth',1)
end

xlabel('Passage Time (seconds)','FontName','Gill Sans','fontsize',14);
ylabel('First-passage probability density','FontName','Gill Sans','fontsize',14);

xlim(axes1,[0.01 10000]);
ylim(axes1,[1E-7 2]);
box(axes1,'on');
set(axes1,'FontName','Gill Sans','FontSize',14,'LineWidth',1.2,'XMinorTick',...
    'on','XScale','log','YMinorTick','on','YScale','log');
set(f5,'Position',[200 300 560 420])
set(gcf, 'Color', 'w');

%% output parameters
str_array_A = {'rf 1 exp','rb 1 exp','k1 1 exp','p1 1 exp','p2 1 exp'};
str_array_B = {'rf 2 exp','rb 2 exp','k1 2 exp','k2 2 exp','p1 2 exp','p2 2 exp','p3 2 exp'};
str_array_B3 = {'rf 3 exp','rb 3 exp','k1 3 exp','k2 3 exp','k3 3 exp','p1 3 exp','p2 3 exp','p3 3 exp','p4 3 exp'};
str_array_B4 = {'rf 4 exp','rb 4 exp','k1 4 exp','k2 4 exp','k3 4 exp','k4 4 exp','p1 4 exp','p2 4 exp','p3 4 exp','p4 4 exp','p5 4 exp'};

[aic,bic] = aicbic([bic_one,bic_two,bic_three,bic_four],[5,7,9,11],[length(pass_times),length(pass_times),length(pass_times),length(pass_times)])

f3 = figure(3);
xlim([0 1]);
ylim([0 1]);
set(gca,'visible','off')
set(gcf, 'Color', 'w');
set(f3, 'Position', [10   197   1400   590] )

for h = 1:length(str_array_A)
    text(0.0,1-((h-1)*0.05),[str_array_A{h} ' = ' num2str(A(h),'%.4f') ' \pm ' num2str(std_param(h),'%.4f')],'fontsize',14);
end

for z = 1:length(str_array_B)
    text(0.3,1-((z-1)*0.05),[str_array_B{z} ' = ' num2str(B(z),'%.4f') ' \pm ' num2str(std_param(z),'%.4f')],'fontsize',14);
end

for r = 1:length(str_array_B3)
    text(0.6,1-((r-1)*0.05),[str_array_B3{r} ' = ' num2str(C(r),'%.4f')],'fontsize',14);
end

for r = 1:length(str_array_B3)
    text(0.8,1-((r-1)*0.05),[str_array_B4{r} ' = ' num2str(D(r),'%.4f')],'fontsize',14);
end

text(0.0,0.2,['aic single = ',num2str(aic(1),'%.2f')],'fontsize',14);
text(0.3,0.2,['aic double = ',num2str(aic(2),'%.2f')],'fontsize',14);
text(0.6,0.2,['aic 3 = ',num2str(aic(3),'%.2f')],'fontsize',14);
text(0.8,0.2,['aic 4 = ',num2str(aic(4),'%.2f')],'fontsize',14);

text(0.0,0.15,['bic single = ',num2str(bic(1),'%.2f')],'fontsize',14);
text(0.3,0.15,['bic double = ',num2str(bic(2),'%.2f')],'fontsize',14);
text(0.6,0.15,['bic 3 = ',num2str(bic(3),'%.2f')],'fontsize',14);
text(0.8,0.15,['bic 4 = ',num2str(bic(4),'%.2f')],'fontsize',14);
end

function out = LLexp_and_LLRW(A,m,tau,interval_bp)
% single exponential + random walk log likelihood function
first = (interval_bp./tau)*((A(1)/A(2))^(interval_bp/2)).*exp(((-((A(1) + A(2))/1))+(( (A(1) + A(2))*sqrt(1-(((A(1) - A(2))/(A(1) + A(2)))^2)) ))).*tau);
second = besseli(interval_bp, ((A(1) + A(2))/1)*sqrt(1-(((A(1) - A(2))/(A(1) + A(2)))^2))*tau,1    );
third = A(3)*exp(-A(3).*tau);

for j = 1:m
    test(j) = (A(4).*(first(j).*second(j))) + (A(5).*third(j));
end

for j = 1:length(test)
    test2(j) = log(test(j));
end

out = sum(test2);

end

function out = LLdoubleexp_and_LLRW(B,m,tau,interval_bp)
% double exponential + random walk log likelihood function
first = (interval_bp./tau)*((B(1)/B(2))^(interval_bp/2)).*exp(((-((B(1) + B(2))/1))+(( (B(1) + B(2))*sqrt(1-(((B(1) - B(2))/(B(1) + B(2)))^2)) ))).*tau);
second = besseli(interval_bp, ((B(1) + B(2))/1)*sqrt(1-(((B(1) - B(2))/(B(1) + B(2)))^2))*tau,1    );
third = B(3)*exp(-B(3).*tau);
fourth = B(4)*exp(-B(4).*tau);

for j = 1:m
    test(j) = (B(5).*(first(j).*second(j))) + (B(6).*third(j)) + (B(7).*fourth(j));
end

for j = 1:length(test)
    test2(j) = log(test(j));
end

out = sum(test2);

end

function out = LL_3_exp_and_LLRW(C,m,tau,interval_bp)
% triple exponential + random walk log likelihood function
first = (interval_bp./tau)*((C(1)/C(2))^(interval_bp/2)).*exp(((-((C(1) + C(2))/1))+(( (C(1) + C(2))*sqrt(1-(((C(1) - C(2))/(C(1) + C(2)))^2)) ))).*tau);
second = besseli(interval_bp, ((C(1) + C(2))/1)*sqrt(1-(((C(1) - C(2))/(C(1) + C(2)))^2))*tau,1    );
third = C(3)*exp(-C(3).*tau);
fourth = C(4)*exp(-C(4).*tau);
fifth = C(5)*exp(-C(5).*tau);

for j = 1:m
    test(j) = (C(6).*(first(j).*second(j))) + (C(7).*third(j)) + (C(8).*fourth(j)) + (C(9).*fifth(j));
end

for j = 1:length(test)
    test2(j) = log(test(j));
end

out = sum(test2);

end

function out = LL_4_exp_and_LLRW(C,m,tau,interval_bp)
% quadrupel exponential + random walk log likelihood function
first = (interval_bp./tau)*((C(1)/C(2))^(interval_bp/2)).*exp(((-((C(1) + C(2))/1))+(( (C(1) + C(2))*sqrt(1-(((C(1) - C(2))/(C(1) + C(2)))^2)) ))).*tau);
second = besseli(interval_bp, ((C(1) + C(2))/1)*sqrt(1-(((C(1) - C(2))/(C(1) + C(2)))^2))*tau,1    );
third = C(3)*exp(-C(3).*tau);
fourth = C(4)*exp(-C(4).*tau);
fifth = C(5)*exp(-C(5).*tau);
sixth = C(6)*exp(-C(6).*tau);

for j = 1:m
    test(j) = (C(7).*(first(j).*second(j))) + (C(8).*third(j)) + (C(9).*fourth(j)) + (C(10).*fifth(j)) + (C(11).*sixth(j));
end

for j = 1:length(test)
    test2(j) = log(test(j));
end

out = sum(test2);

end