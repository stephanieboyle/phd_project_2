%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Frequency Analysis - for individual subjects
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subs = dir('W:\Project 1\log\EEG Data\Subj_*');
addpath('W:\Project 1\log\EEG Data\functions'); 
cd('W:\Project 1\log\EEG Data\Analysis_TFR')
addpath('W:\Project 1\log\EEG Data\functions')

% Loop through subjects
for subj = 1:length(subs);
    
    cd(sprintf('W:/Project 1/log/EEG Data/%s',subs(subj).name))
    files = dir('PreproICA_S*_EOG_out.mat');
    TFR = cell(1,length(files));
    
    for k = 1:length(files)
        
        load(files(k).name)
        twn = (-0.5:0.05:1.2);
        baselinetype = 'relative';
        
        % lower frequencies config
        cfg_tfmt              = [];
        cfg_tfmt.toi          = twn;
        cfg_tfmt.trials = 'all';
        cfg_tfmt.output       = 'pow';
        cfg_tfmt.channel      = 'all';
        cfg_tfmt.method       = 'mtmconvol';
        cfg_tfmt.foi          = [1:10,12:2:24];
        cfg_tfmt.taper        = 'hanning';
        cfg_tfmt.keeptrials = 'yes';
        cfg_tfmt.t_ftimwin    = 5./cfg_tfmt.foi;            % frequency depending window
        tmp = ft_freqanalysis(cfg_tfmt, dataX);
        
        % multiple tapers above 24Hz config
        cfg_tfmt2              = [];
        cfg_tfmt2.trials = 'all';
        cfg_tfmt2.toi          = twn;
        cfg_tfmt2.output       = 'pow';
        cfg_tfmt2.channel      = 'all';
        cfg_tfmt2.method       = 'mtmconvol';
        cfg_tfmt2.foi          = (26:4:80);
        cfg_tfmt2.t_ftimwin    = 5./cfg_tfmt2.foi;
        cfg_tfmt2.keeptrials = 'yes';
        cfg_tfmt2.tapsmofrq  = 0.4 *cfg_tfmt2.foi;          % frequency smoothing
        
        % high frequencies
        tmp2 = ft_freqanalysis(cfg_tfmt2, dataX);
        
        % append both frequency parts together
        cfg =[];
        cfg.appenddim = 'freq';
        cfg.parameter  = 'powspctrm';
        combined = ft_appendfreq(cfg, tmp, tmp2);
        combined.trialinfo = dataX.trialinfo ; % save trialinfo;
        
        % normalize to baseline
        block = 1;
        cfg =[];
        cfg.baseline     = [-0.5 -0.1];
        cfg.baselinetype = baselinetype;
        TFR{k} = ft_freqbaseline(cfg, combined);
        
    end
    
    % put all TFR blocks together into one TFR structure
    cfg = [];
    cfg.parameter = 'powspctrm';
    TFR_all = ft_appendfreq(cfg,TFR{:});
    
    % save the file
    sname = sprintf('S0%d_TFR_EOG',subj);
    save(sname,'TFR_all')
    fprintf('TFR S0%d \n',subj)
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Group analysis : AVH and AVL separately 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get AVH-VH, AVL-VL powspctrms
clear; clc; close all
cd('W:/Project 1/log/EEG Data/Analysis_TFR');
addpath('W:\Project 1\log\EEG Data\functions');
subs = dir('Subj_*');
tfHigh = zeros(20,64,31,35);                            % AVH
tfLow = zeros(20,64,31,35);                             % AVL 

% Loop through subjects
for subj=1:length(subs);
    
    cd(sprintf('W:/Project 1/log/EEG Data/%s',subs(subj).name))
    load(sprintf('S0%d_TFR_EOG.mat',subj))
    
    if subj==1;
        tax = TFR_all.time; fax = TFR_all.freq;
    end
    
    % get conditions:
    Df = cell(1,9);
    for k = 1:9;
        Df{k} = find(TFR_all.trialinfo(:,3)==k);
    end
    
    % pwspctrms:
    tfHigh(subj,:,:,:) = sq(mean(TFR_all.powspctrm(Df{4},1:64,:,:))) - sq(mean(TFR_all.powspctrm(Df{1},1:64,:,:))); % AVH-VH
    tfLow(subj,:,:,:) = sq(mean(TFR_all.powspctrm(Df{7},1:64,:,:))) - sq(mean(TFR_all.powspctrm(Df{2},1:64,:,:)));  % AVL-VL

    fprintf('S0%d done \n',subj);
    clear TFR_all
end

save('TF_High.mat','tfHigh','tax','fax')
save('TF_Low.mat','tfLow','tax','fax')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time frequency for each condition separately
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all
cd('W:/Project 1/log/EEG Data');
subs = dir('Subj_*');
VH = zeros(67,31,35,length(subs));
VL = VH; AUD = VH; AVH = VH; AVL = VH;              % 0 conflict
AVH2 = VH; AVH_2 = VH; AVL2 = VH; AVL_2 = VH;       % +- conflict 

for s = 1:length(subs);
    
    cd('W:\Project 1\log\EEG Data\Analysis_TFR\Individual Subs')
    load(sprintf('S0%d_TFR_EOG.mat',s))
    
    % get conditions:
    Df = cell(1,9);
    for k = 1:9;
        Df{k} = find(TFR_all.trialinfo(:,3)==k);
    end
    
    % Unisensory conditions 
    VH(:,:,:,s) = sq(mean(TFR_all.powspctrm(Df{1},:,:,:)));         % visual high
    VL(:,:,:,s) = sq(mean(TFR_all.powspctrm(Df{2},:,:,:)));         % visual low 
    AUD(:,:,:,s) = sq(mean(TFR_all.powspctrm(Df{3},:,:,:)));        % auditory 
    
    % No conflict conditions
    AVH(:,:,:,s) = sq(mean(TFR_all.powspctrm(Df{4},:,:,:)));        % audiovisual high
    AVL(:,:,:,s) = sq(mean(TFR_all.powspctrm(Df{7},:,:,:)));        % audiovisual low
    
    % Conflict conditions 
    AVH2(:,:,:,s) = sq(mean(TFR_all.powspctrm(Df{5},:,:,:)));       % audiovisual high +2
    AVH_2(:,:,:,s) = sq(mean(TFR_all.powspctrm(Df{6},:,:,:)));      % audiovisual high -2
    AVL2(:,:,:,s) = sq(mean(TFR_all.powspctrm(Df{8},:,:,:)));       % audiovisual low +2
    AVL_2(:,:,:,s) = sq(mean(TFR_all.powspctrm(Df{9},:,:,:)));      % audiovisual low -2

    fprintf('S0%d done... \n',s)
end
cd('W:\Project 1\log\EEG Data\Analysis_TFR'); 
load('Group_TFR','fax','tax')


% save stuff 
save('Group_AUD.mat','AUD','fax','tax'); 
save('Group_VH.mat','VH','fax','tax'); 
save('Group_VL.mat','VL','fax','tax'); 
save('Group_AVH.mat','AVH','fax','tax'); 
save('Group_AVL.mat','AVL','fax','tax'); 
save('Group_AVH2.mat','AVH2','fax','tax');
save('Group_AVH_2.mat','AVH_2','fax','tax');
save('Group_AVL2.mat','AVL2','fax','tax');
save('Group_AVL_2.mat','AVL_2','fax','tax');

% Average over the group and electrodes
Group(:,:,:,1) = sq(mean(VH(1:64,:,:,:),4));                % visual high
Group(:,:,:,2) = sq(mean(VL(1:64,:,:,:),4));                % visual low 
Group(:,:,:,3) = sq(mean(AUD(1:64,:,:,:),4));               % auditory 
Group(:,:,:,4) = sq(mean(AVH(1:64,:,:,:),4));               % audiovisual high
Group(:,:,:,5) = sq(mean(AVL(1:64,:,:,:),4));               % audiovisual low 

GroupConf(:,:,:,1) = sq(mean(AVH2(1:64,:,:,:),4));          % audiovisual high +2
GroupConf(:,:,:,2) = sq(mean(AVH_2(1:64,:,:,:),4));         % audiovisual high -2
GroupConf(:,:,:,3) = sq(mean(AVL2(1:64,:,:,:),4));          % audiovisual low +2
GroupConf(:,:,:,4) = sq(mean(AVL_2(1:64,:,:,:),4));         % audiovisual low -2 

% save stuff 
names = {'VH','VL','AUD','AVH','AVL'};
namesC = {'AVH2','AVH_2','AVL','AVL_2'}; 
save('GROUP_TFR.mat','Group','GroupConf','tax','fax','names','namesC')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Group TFR representations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
clim = [-0.5 5]; 
tR = find(fax<50); 

GroupAvg = sq(mean(Group,1));
GroupConfAvg = sq(mean(GroupConf,1)); 

% Plot unisensory trials 
for k = 1:3; 
    subplot(3,3,k);
    imagesc(tax,fax(tR),sq(GroupAvg(tR,:,k)),clim); colorbar
    title(names{k})
end

% plot AVH trials 
subplot 334
imagesc(tax,fax(tR),sq(GroupAvg(tR,:,4)),clim); colorbar; title('AVH');
subplot 335;
imagesc(tax,fax(tR),sq(GroupConfAvg(tR,:,1)),clim); colorbar; title('AVH +2')
subplot 336;
imagesc(tax,fax(tR),sq(GroupConfAvg(tR,:,2)),clim); colorbar; title('AVH -2')

% plot AVL trials 
subplot 337
imagesc(tax,fax(tR),sq(GroupAvg(tR,:,5)),clim); colorbar; title('AVL');
subplot 338;
imagesc(tax,fax(tR),sq(GroupConfAvg(tR,:,3)),clim); colorbar; title('AVL +2')
subplot 339;
imagesc(tax,fax(tR),sq(GroupConfAvg(tR,:,4)),clim); colorbar; title('AVL -2')

suptitle('Group Averaged over Subs and Electrodes')
saveas(gcf,'Figure_GroupTFR.fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Work out the double differences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AVH_VH = Group(:,tR,:,4) - Group(:,tR,:,1);         % AVH-VH (frequencies below 50 Hz);
AVL_VL = Group(:,tR,:,5) - Group(:,tR,:,2);         % AVL-VL (frequencies below 50 Hz); 
AV_DD = AVL_VL - AVH_VH;                            % [AVL-VL] - [AVH-VH]; 

clim = [-0.25 0.5];

% plot high reliability difference
subplot 131
imagesc(tax,fax(tR),sq(mean(mean(AVH_VH,4),1)),clim); colorbar; 
title('AVH-VH'); 

% plot low reliability difference
subplot 132
imagesc(tax,fax(tR),sq(mean(mean(AVL_VL,4),1)),clim); colorbar; 
title('AVL-VL'); 

% plot low-high reliability difference
subplot 133
imagesc(tax,fax(tR),sq(mean(mean(AV_DD,4),1)),clim); colorbar; 
title('(AVL-VL) - (AVH-VH)'); 

suptitle('Group TFR Differences')
saveas(gcf,'Figure_GroupTFR_Differences.fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multiplots for the differences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('TF64'); 
TF64.freq = TF64.freq(tR);                          % just take the narrower range
TF_AVH_VH = TF64;
TF_AVH_VH.powspctrm = AVH_VH;                       % AVH-VH structure
TF_AVL_VL = TF64;
TF_AVL_VL.powspctrm = AVL_VL;                       % AVL-VL structure

% Make the multiplots 
cfg = []; 
cfg.layout = 'biosemi64.lay'; 
cfg.zlim = [-0.25 0.45]; 
subplot 121
ft_multiplotTFR(cfg, TF_AVH_VH); colorbar           % AVH-VH
title('AVH-VH'); 
subplot 122
ft_multiplotTFR(cfg, TF_AVL_VL); colorbar           % AVL-VL
title('AVL-VL'); 

saveas(gcf,'Figure_GroupTFR_DifferencesMTP.fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do imagesc for each of the four electrode groups to see where the
% differences are more clearly 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[electrodes,Index,Ylab] = eegsb_GroupLabels(); 
clf
clim = [-0.25 0.25; -1.5 0.8 ; -0.2 0.6 ; -0.35 0.5]; 
for k = 1:4; 
    
    % plot the different electrode groupings for AVH-VH 
    subplot(3,4,k)
    imagesc(tax,fax(tR),sq(mean(AVH_VH(Index{k},:,:))),clim(k,:)); colorbar
    title(Ylab{k}); ylabel('Freq Hz'); xlabel('Time(s)');
    np = np+1; 
    if k==1; ylabel(sprintf('AVH - VH TRIALS \n Freq Hz'),'FontSize',10); end
    
    % plot the different electrode groupings for AVL-VL 
    subplot(3,4,k+4)
    imagesc(tax,fax(tR),sq(mean(AVL_VL(Index{k},:,:))),clim(k,:)); colorbar
    title(Ylab{k}); ylabel('Freq Hz'); xlabel('Time(s)');
    if k==1; ylabel(sprintf('AVL - VL TRIALS \n Freq Hz'),'FontSize',10); end
    
    % Plot the double difference 
    subplot(3,4,k+8)
    imagesc(tax,fax(tR),sq(mean(AV_DD(Index{k},:,:))),clim(k,:)); colorbar
    title(Ylab{k}); ylabel('Freq Hz'); xlabel('Time(s)');
    if k==1; ylabel(sprintf('[AVL-VL]-[AVH-VH] TRIALS \n Freq Hz'),'FontSize',10); end
        
end

saveas(gcf,'Figure_GroupTFR_Differences_ElectrodeGroup.fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% From the plots, it seems that there are differences in low power bands
% (theta) and differences in beta (12-15) early in the trials, and perhaps
% in the higher bands (>30Hz) later in the trial. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster analysis : [AVL-VH] vs [AVH-VH] THETA 4-6 HZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[electrodes,Index,Ylab] = eegsb_GroupLabels(); 
load('Group_AVH.mat')
load('Group_AVL.mat'); 
load('Group_VH.mat'); 
load('Group_VL.mat'); 
load('TF64')

% AVH
gAVH = TF64; 
gAVH.dimord = 'subj_chan_freq_time';
gAVH.powspctrm = permute(AVH(1:64,:,:,:),[4,1,2,3]);

% AVL 
gAVL = gAVH; 
gAVL.powspctrm = permute(AVL(1:64,:,:,:),[4,1,2,3]);

% AVH-VH
group_AVH = TF64; 
group_AVH.dimord = 'subj_chan_freq_time';
group_AVH.powspctrm = AVH(1:64,:,:,:) - VH(1:64,:,:,:);         % AVH-VH
group_AVH.powspctrm = permute(group_AVH.powspctrm,[4,1,2,3]); 

% AVL-VL
group_AVL = group_AVH; 
group_AVL.powspctrm = AVL(1:64,:,:,:) - VL(1:64,:,:,:);         % AVL-VL 
group_AVL.powspctrm = permute(group_AVL.powspctrm,[4,1,2,3]); 

cfg_neighb = eegck_BiosemiNbrs(TF64);                           % Define neighbours 
nsub = size(group_AVH.powspctrm,1);                             % number of subs 

% Set up the different things you want: we know that in the main
% experiment, our effects were all happening before 300ms, so we will
% restrict it to that for now and won't average over time. 
cfg = [];
cfg.channel = 'all';
cfg.frequency        = [4 6];                   % lower frequency band (theta)
cfg.latency          = [0 0.3];                 % keep the early time window. 
cfg.parameter = 'powspctrm'; 
cfg.avgoverfreq = 'yes';
cfg.statistic        = 'ft_statfun_depsamplesT'; % t-test
cfg.tail             = 0;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';      % CORRECT FOR 2 TAILED !!!!!!
cfg.neighbours       = cfg_neighb;
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.025;
cfg.clusterstatistic = 'maxsize';
cfg.clustertail      = 0;
cfg.minnbchan        = 2;
cfg.numrandomization = 1000;

% design matrix;
design = zeros(2,2*nsub);
design(1,:) = [(1:nsub),(1:nsub)]; % identifies subjects (uvar)
design(2,:) = [ones(1,nsub)*1,ones(1,nsub)*2]; % identifies conditions (ivar)
cfg.design   = design;
cfg.uvar     = 1; % unit of observation (subjects)
cfg.ivar     = 2; % condition
cfg.cvar =[];

% time -frequency graph to compare conditions
[stat_theta] = ft_freqstatistics(cfg,group_AVL,group_AVH);

% collect cluster informatiocn
[neg_signif_clust,pos_signif_clust,neg,pos,stat_neg,stat_pos] = eegck_get_clusters(stat_theta);

% Plot where the significant values are 
figure
imagesc(cfg.latency,1:64,sq(neg)); colorbar
title(sprintf('Freq %s:%s Hz, AVL Diff vs. AVH Diff \n Avdg over freq',num2str(cfg.frequency(1)),num2str(cfg.frequency(2)))); 

% Plot the tvalues in a topoplot
load('FtDummy','EvpDummy');
EvpDummy.time = stat_theta.time; 
EvpDummy.avg = stat_theta.stat; 

clf; 
cfg_topo = []; 
cfg_topo.highlight = 'on';
cfg_topo.zlim = [-3 3];
cfg_topo.highlightsymbol = '*';
cfg_topo.highlightcolor = [1 1 1];
cfg_topo.highlightsize = 5;
cfg_topo.comment = ' ';
cfg_topo.layout = 'biosemi64.lay'; 
% cfg_topo.xlim = stat_theta.time; 
% cfg_topo.highlightchannel = find(neg(:,1,:)~=0); 
% EvpDummy.avg = stat_theta.stat; 

for k = 1:size(stat_theta.stat,3);
    clf
    EvpDummy.time = stat_theta.time(k); 
    EvpDummy.avg =  sq(stat_theta.stat(:,:,k));
    cfg_topo.highlightchannel    = find(neg(:,1,k)~=0);
    ft_topoplotER(cfg_topo, EvpDummy);
    colorbar
    suptitle(sprintf('Theta T values Time: %s',num2str(stat_theta.time(k))));
    colormap(redblue)
    sname = sprintf('Figure_TVALS_Theta_%sms.fig',num2str(stat_theta.time(k)));
    saveas(gcf,sname)
%     suptitle('THETA T VALUES'); colormap(redblue); 
end




% plot the raw effect
fax = group_AVL.freq;                                                       % all frequencies
tax = group_AVL.time;                                                       % all time points 
fUse = find(fax>=cfg.frequency(1) & fax<=cfg.frequency(end));               % frequencies used
tUse = find(tax>=cfg.latency(1) & tax<=cfg.latency(end));                   % times used
raweffectAVL = sq(mean(mean(group_AVL.powspctrm(:,:,fUse,tUse),1),3));      % AVL raw effect
raweffectAVH = sq(mean(mean(group_AVH.powspctrm(:,:,fUse,tUse),1),3));      % AVH raw effect
rawDiff = raweffectAVL - raweffectAVH;

% plot raw effect
clf
cfg_topo = []; 
cfg_topo.highlight = 'on';
cfg_topo.zlim = [-0.2 0.2];
cfg_topo.highlightsymbol = '*';
cfg_topo.highlightcolor = [1 1 1];
cfg_topo.highlightsize = 5;
cfg_topo.comment = ' ';
cfg_topo.layout = 'biosemi64.lay'; 

for k = 1:size(stat_theta.stat,3);
    clf
    EvpDummy.avg =  sq(rawDiff(:,k));
    cfg_topo.highlightchannel    = find(neg(:,1,k)~=0);
    ft_topoplotER(cfg_topo, EvpDummy);
    colorbar
    title(sprintf('Theta Raw vals [AVL-VL]-[AVH-VH], Time: %s',num2str(stat_theta.time(k))));
    sname = sprintf('Figure_RAW_Theta_%sms.fig',num2str(stat_theta.time(k)));
    saveas(gcf,sname)
end
% suptitle('[AVL-VL] - [AVH-VH] raw diff VALUES'); colormap(jet); 
% saveas(gcf,'Figure_RAW_Theta.fig')

save('Stats_ThetaHz.mat','stat_theta','group_AVL','group_AVH','cfg','cfg_topo','rawDiff')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
% average the significnat time window 
load('FtDummy','EvpDummy');
EvpDummy.time = stat_theta.time(1:4); 
EvpDummy.avg = sq(stat_theta.stat(:,1,1:4)); 

clf; 
cfg_topo = []; 
cfg_topo.highlight = 'on';
cfg_topo.zlim = [-3 3];
cfg_topo.highlightsymbol = '*';
cfg_topo.highlightcolor = [1 1 1];
cfg_topo.highlightsize = 5;
cfg_topo.comment = ' ';
cfg_topo.layout = 'biosemi64.lay'; 
cfg_topo.highlightchannel = find(sum(sq(neg),2)~=0); 
ft_topoplotER(cfg_topo, EvpDummy);
colormap(redblue); colorbar

% plot the raw effect
fax = group_AVL.freq;                                                       % all frequencies
tax = group_AVL.time;                                                       % all time points 
fUse = find(fax>=cfg.frequency(1) & fax<=cfg.frequency(end));               % frequencies used
tUse = find(tax>=cfg.latency(1) & tax<=cfg.latency(end));                   % times used
raweffectAVL = sq(mean(mean(group_AVL.powspctrm(:,:,fUse,tUse),1),3));      % AVL raw effect
raweffectAVH = sq(mean(mean(group_AVH.powspctrm(:,:,fUse,tUse),1),3));      % AVH raw effect
rawDiff = raweffectAVL - raweffectAVH;

% plot raw effect
clf
cfg_topo = []; 
cfg_topo.highlight = 'on';
cfg_topo.zlim = [-0.2 0.2];
cfg_topo.highlightsymbol = '*';
cfg_topo.highlightcolor = [1 1 1];
cfg_topo.highlightsize = 5;
cfg_topo.comment = ' ';
cfg_topo.layout = 'biosemi64.lay'; 
EvpDummy.avg =  sq(rawDiff(:,1:4));
cfg_topo.highlightchannel    = find(sum(sq(neg),2)~=0); 
ft_topoplotER(cfg_topo, EvpDummy);
colorbar; colormap(jet)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cluster analysis : [AVL-VH] vs [AVH-VH] BETA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg_neighb = eegck_BiosemiNbrs(TF64);                           % Define neighbours 
nsub = size(group_AVH.powspctrm,1);                             % number of subs 

% Set up the different things you want: we know that in the main
% experiment, our effects were all happening before 300ms, so we will
% restrict it to that for now and won't average over time. 
cfg = [];
cfg.channel = 'all';
cfg.frequency        = [10 12];
cfg.latency          = [0.1 0.4];
cfg.parameter = 'powspctrm'; 
cfg.avgoverfreq = 'yes';
cfg.statistic        = 'ft_statfun_depsamplesT'; % t-test
cfg.tail             = 0;
cfg.alpha            = 0.05;
cfg.correcttail      = 'alpha';      % CORRECT FOR 2 TAILED !!!!!!
cfg.neighbours       = cfg_neighb;
cfg.method           = 'montecarlo';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.025;
cfg.clusterstatistic = 'maxsize';
cfg.clustertail      = 0;
cfg.minnbchan        = 2;
cfg.numrandomization = 1000;

% design matrix;
design = zeros(2,2*nsub);
design(1,:) = [(1:nsub),(1:nsub)]; % identifies subjects (uvar)
design(2,:) = [ones(1,nsub)*1,ones(1,nsub)*2]; % identifies conditions (ivar)
cfg.design   = design;
cfg.uvar     = 1; % unit of observation (subjects)
cfg.ivar     = 2; % condition
cfg.cvar =[];

% time -frequency graph to compare conditions
[stat_beta] = ft_freqstatistics(cfg,group_AVL,group_AVH);

% collect cluster informatiocn
[neg_signif_clust,pos_signif_clust,neg,pos,stat_neg,stat_pos] = eegck_get_clusters(stat_beta);

% Plot where the significant values are 
figure
imagesc(cfg.latency,1:64,sq(pos)); colorbar
title(sprintf('Freq %s:%s Hz, AVL Diff vs. AVH Diff \n Avdg over freq',num2str(cfg.frequency(1)),num2str(cfg.frequency(2)))); 

% Plot the tvalues in a topoplot
load('FtDummy','EvpDummy');
EvpDummy.time = stat_beta.time; 
EvpDummy.avg = stat_beta.stat; 

clf; 
cfg_topo = []; 
cfg_topo.highlight = 'on';
cfg_topo.zlim = [-3 3];
cfg_topo.highlightsymbol = '*';
cfg_topo.highlightcolor = [1 1 1];
cfg_topo.highlightsize = 5;
cfg_topo.comment = ' ';
cfg_topo.layout = 'biosemi64.lay'; 

for k = 5; %1:size(stat_beta.stat,3);
%     subplot(3,3,k); hold on
    EvpDummy.time = stat_beta.time(k); 
    EvpDummy.avg =  sq(stat_beta.stat(:,:,k));
    cfg_topo.highlightchannel    = find(pos(:,1,k)~=0);
    ft_topoplotER(cfg_topo, EvpDummy);
    colorbar
    title(sprintf('Time: %s',num2str(stat_beta.time(k))));
end
suptitle('BETA T VALUES'); colormap(redblue); 
saveas(gcf,'Figure_TVALS_Beta.fig')


% plot the raw effect
fax = group_AVL.freq;                                                       % all frequencies
tax = group_AVL.time;                                                       % all time points 
fUse = find(fax>=cfg.frequency(1) & fax<=cfg.frequency(end));               % frequencies used
tUse = find(tax>=cfg.latency(1) & tax<=cfg.latency(end));                   % times used
raweffectAVL = sq(mean(mean(group_AVL.powspctrm(:,:,fUse,tUse),1),3));      % AVL raw effect
raweffectAVH = sq(mean(mean(group_AVH.powspctrm(:,:,fUse,tUse),1),3));      % AVH raw effect
rawDiff = raweffectAVL - raweffectAVH;

% plot raw effect
clf
cfg_topo = []; 
cfg_topo.highlight = 'on';
cfg_topo.zlim = [-0.2 0.2];
cfg_topo.highlightsymbol = '*';
cfg_topo.highlightcolor = [1 1 1];
cfg_topo.highlightsize = 5;
cfg_topo.comment = ' ';
cfg_topo.layout = 'biosemi64.lay'; 

for k = 1:size(stat_beta.stat,3);
    subplot(3,3,k); hold on
    EvpDummy.avg =  sq(rawDiff(:,k));
    cfg_topo.highlightchannel    = find(pos(:,1,k)~=0);
    ft_topoplotER(cfg_topo, EvpDummy);
    colorbar
    title(sprintf('Time: %s',num2str(stat_beta.time(k))));
end
suptitle('[AVL-VL] - [AVH-VH] raw diff BETA VALUES'); colormap(jet); 
saveas(gcf,'Figure_RAW_Beta.fig')

save('Stats_BetaHz.mat','stat_beta','group_AVL','group_AVH','cfg','cfg_topo','rawDiff')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% results : theta significance and beta early after stimulus onset 
% theta first, beta second. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Try Gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NO GAMMA




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Look at individual conditions power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [electrodes,Index,Ylab] = eegsb_GroupLabels(); 
load('Group_AVH.mat')
load('Group_AVL.mat'); 
load('Group_VH.mat'); 
load('Group_VL.mat'); 
load('Group_AUD.mat'); 

% Work out the diffs
AVL_VL = AVL-VL;
AVH_VH = AVH-VH; 

% initialise things 
Rows = eegck_rowlabels;
e = 1:64; 
theta = find(fax>=4 & fax<=6); 
beta = find(fax>=10 & fax<=12); 
t = find(tax>=0 & tax<=0.3); 
climT =[1.5 2.7];
climB = [1 2];

%--------------------------------------------------------------------------
% Unisensory plots
%--------------------------------------------------------------------------

% Auditory 
subplot 321
imagesc(tax(t),1:64,sq(mean(mean(AUD(Rows.RowLabels2,theta,t,:),4),2)),climT)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName);
title('theta 4:6 Hz'); ylabel('Auditory','FontSize',12); colorbar
subplot 322
imagesc(tax(t),1:64,sq(mean(mean(AUD(Rows.RowLabels2,beta,t,:),4),2)),climB)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName)
title('beta 10:12 Hz'); colorbar

% Visual High
subplot 323
imagesc(tax(t),1:64,sq(mean(mean(VH(Rows.RowLabels2,theta,t,:),4),2)),climT)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName);
title('theta 4:6 Hz'); ylabel('Visual High','FontSize',12); colorbar
subplot 324
imagesc(tax(t),1:64,sq(mean(mean(VH(Rows.RowLabels2,beta,t,:),4),2)),climB)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName)
title('beta 10:12 Hz'); colorbar

% Visual Low 
subplot 325
imagesc(tax(t),1:64,sq(mean(mean(VL(Rows.RowLabels2,theta,t,:),4),2)),climT)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName);
title('theta 4:6 Hz'); ylabel('Visual Low','FontSize',12); colorbar
subplot 326
imagesc(tax(t),1:64,sq(mean(mean(VL(Rows.RowLabels2,beta,t,:),4),2)),climB)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName)
title('beta 10:12 Hz'); colorbar

suptitle('Unisensory Theta & Beta, Group Averaged Power')
saveas(gcf,'Figure_BetaTheta_Unisensory.fig')


%--------------------------------------------------------------------------
% Mutlisensory plots : audiovisual 
%--------------------------------------------------------------------------
climT =[1 3.5];
climB = [1 2];
clf 

% Audiovisual high  
subplot 221
imagesc(tax(t),1:64,sq(mean(mean(AVH(Rows.RowLabels2,theta,t,:),4),2)),climT)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName);
title('theta'); ylabel('Audiovisual High','FontSize',12); colorbar
subplot 222
imagesc(tax(t),1:64,sq(mean(mean(AVH(Rows.RowLabels2,beta,t,:),4),2)),climB)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName)
title('beta'); colorbar

% Audiovisual low
subplot 223
imagesc(tax(t),1:64,sq(mean(mean(AVL(Rows.RowLabels2,theta,t,:),4),2)),climT)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName);
title('theta'); ylabel('Audiovisual Low','FontSize',12); colorbar
subplot 224
imagesc(tax(t),1:64,sq(mean(mean(AVL(Rows.RowLabels2,beta,t,:),4),2)),climB)
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName)
title('beta'); colorbar

suptitle('Audiovisual Theta & Beta, Group Averaged Power')
saveas(gcf,'Figure_BetaTheta_Audiovisual.fig')


%--------------------------------------------------------------------------
% Mutlisensory plots : audiovisual difference
%--------------------------------------------------------------------------
AV_Diff = AVL_VL-AVH_VH;

% figure
clf; subplot 121
imagesc(tax(t),1:64,sq(mean(mean(AV_Diff(Rows.RowLabels2,theta,t,:),4),2)))
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName);
title('theta'); ylabel('Audiovisual High','FontSize',12); colorbar
subplot 122
imagesc(tax(t),1:64,sq(mean(mean(AV_Diff(Rows.RowLabels2,beta,t,:),4),2)))
set(gca,'YTick',[1 cumsum(Rows.RowL(1:end-1))+1],'YTickLabel',Rows.RowName)
title('beta'); colorbar

suptitle('[AVL-VL] - [AVH-VH], Theta & Beta, Group Averaged Power')
saveas(gcf,'Figure_BetaTheta_AVdifference.fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Power line graphs for each subject and group mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = find(tax>=0 & tax<=0.5);                        % shorter time window
betaPower = sq(mean(AV_Diff(e,beta,t,:),2));        % beta power
bPower = mean(betaPower,3);                         % group averaged beta power
thetaPower = sq(mean(AV_Diff(e,theta,t,:),2));      % theta power
tPower = mean(thetaPower,3);                        % group averaged theta power

clf
% Theta Power
subplot 121
plot(tax(t),sq(mean(thetaPower)),'color',[1 0.8 0.8])
hold on
plot(tax(t),mean(tPower),'r','LineWidth',3); ylim([-0.6 0.6])
hline(0,':k'); title('Theta Power'); xlim([0 0.3])

% Beta Power
subplot 122
plot(tax(t),sq(mean(betaPower)),'color',[0.8 0.9 1])
hold on
plot(tax(t),mean(bPower),'b','LineWidth',3); ylim([-0.6 0.6])
hline(0,':k'); title('Beta Power'); xlim([0 0.3])

saveas(gcf,'Figure_BetaTheta_PowerPlots_Short.fig')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the theta and beta for AVL-VL and AVH-VH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AVLVL_beta = sq(AVL_VL(1:64,beta,t,:));
AVLVL_theta = sq(AVL_VL(1:64,theta,t,:)); 

AVHVH_beta = sq(AVH_VH(1:64,beta,t,:)); 
AVHVH_theta = sq(AVH_VH(1:64,theta,t,:)); 

betaVals = [sq(mean(mean(AVLVL_beta,1),2)) sq(mean(mean(AVHVH_beta,1),2))]; 
thetaVals = [sq(mean(mean(mean(AVLVL_theta,1),2),3)) sq(mean(mean(mean(AVHVH_theta,1),2),3))]; 

clf; hold on; 
subplot 122
stagger = 0.95:0.005:1.05;
for k =1:20; 
    plot(stagger(k),betaVals(k,1),'o','Color','b'); hold on
    plot(stagger(k)+1,betaVals(k,2),'o','Color','r'); hold on
end
xlim([0 3]); hold on; boxplot(betaVals,'colors','k');
ylabel('Raw Difference Power Vals','FontSize',14); 
set(gca,'XTick',1:2,'XTickLabel',{'AVL-VL','AVH-VH'})
title('Beta frequency band, 0 to 300ms avg','FontSize',14);
ylim([-0.35 0.35]); hline(0,':k')

subplot 121
stagger = 0.95:0.005:1.05;
for k =1:20; 
    plot(stagger(k),thetaVals(k,1),'o','Color','b'); hold on
    plot(stagger(k)+1,thetaVals(k,2),'o','Color','r'); hold on
end
xlim([0 3]); hold on; boxplot(thetaVals,'colors','k');
ylabel('Raw Difference Power Vals','FontSize',14); 
set(gca,'XTick',1:2,'XTickLabel',{'AVL-VL','AVH-VH'})
title('Theta frequency band, 0 to 300ms avg','FontSize',14); 
ylim([-0.35 1.5]);hline(0,':k');
set(gcf,'color','w')
saveas(gcf,'Figure_BetaTheta_Boxplots.fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot the theta and beta for AVL-VL and AVH-VH for the significant points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaT = find(tax>=0 & tax<=0.15);                        % shorter time window
betaT = find(tax==0.3); 

AVLVL_beta = sq(AVL_VL(1:64,beta,betaT,:));
AVLVL_theta = sq(AVL_VL(1:64,theta,thetaT,:)); 

AVHVH_beta = sq(AVH_VH(1:64,beta,betaT,:)); 
AVHVH_theta = sq(AVH_VH(1:64,theta,thetaT,:)); 

betaVals = [sq(mean(AVLVL_beta,1))'  sq(mean(AVHVH_beta,1))'];
thetaVals = [sq(mean(mean(mean(AVLVL_theta,1),2),3)) sq(mean(mean(mean(AVHVH_theta,1),2),3))]; 

% thetaVals = [sq(mean(mean(mean(AVLVL_theta,1),2),3)) sq(mean(mean(mean((AVHVH_theta,1),2),3)))]; 


clf; hold on; 

% THETA
subplot 121
stagger = 0.95:0.005:1.05;
for k =1:20; 
    plot(stagger(k),thetaVals(k,1),'o','Color','b'); hold on
    plot(stagger(k)+1,thetaVals(k,2),'o','Color','r'); hold on
end
xlim([0 3]); hold on; boxplot(thetaVals,'colors','k');
ylabel('Raw Difference Power Vals','FontSize',14); 
set(gca,'XTick',1:2,'XTickLabel',{'AVL-VL','AVH-VH'})
title('Theta frequency band 0 to 150ms avg','FontSize',14); 
ylim([-0.35 1.5]);hline(0,':k');

% BETA
subplot 122
stagger = 0.95:0.005:1.05;
for k =1:20; 
    plot(stagger(k),betaVals(k,1),'o','Color','b'); hold on
    plot(stagger(k)+1,betaVals(k,2),'o','Color','r'); hold on
end
xlim([0 3]); hold on; boxplot(betaVals,'colors','k');
ylabel('Raw Difference Power Vals','FontSize',14); 
set(gca,'XTick',1:2,'XTickLabel',{'AVL-VL','AVH-VH'})
title('Beta frequency band, 300s','FontSize',14);
ylim([-0.4 0.45]); hline(0,':k')

 set(gcf,'color','w')
saveas(gcf,'Figure_BetaTheta_SignificantTP_Boxplots.fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the mean over time for the av conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaE = [6 7 12 13 14 18 19 60]; 
betaE = [7 8 9 12 13 14 18
    
t = find(tax>=0 & tax<=0.5);                               % shorter time window
betaPower_AVLVL = sq(mean(AVL_VL(betaE,beta,t,:),2));             % beta power
bPower_AVLVL = sq(mean(betaPower_AVLVL,1));                         % group averaged beta power
thetaPower_AVLVL = sq(mean(AVL_VL(thetaE,theta,t,:),2));              % theta power
tPower_AVLVL = sq(mean(thetaPower_AVLVL,1));                        % group averaged theta power

betaPower_AVHVH = sq(mean(AVH_VH(betaE,beta,t,:),2));                % beta power
bPower_AVHVH = sq(mean(betaPower_AVHVH,1));                         % group averaged beta power
thetaPower_AVHVH = sq(mean(AVH_VH(thetaE,theta,t,:),2));              % theta power
tPower_AVHVH = sq(mean(thetaPower_AVHVH,1));                        % group averaged theta power

clf
subplot 121
plot(tax(t),tPower_AVLVL,'color',[1 0.85 0.85]); hold on; 
plot(tax(t),tPower_AVHVH,'color',[0.8 0.9 1]); hold on; 
H1 = plot(tax(t),mean(tPower_AVLVL,2),'r','LineWidth',3); hold on
H2 = plot(tax(t),mean(tPower_AVHVH,2),'b','LineWidth',3); 
legend([H1,H2],{'AVL-VL','AVH-VH'})
xlim([0 0.5]); xlabel('Time')
title('Theta Power','FontSize',14); ylabel('Power Vals','FontSize',14); 
ylim([-2.5 2.5]); hline(0,':k')


subplot 122
plot(tax(t),bPower_AVLVL,'color',[1 0.85 0.85]); hold on; 
plot(tax(t),bPower_AVHVH,'color',[0.8 0.9 1]); hold on; 
H1 = plot(tax(t),mean(bPower_AVLVL,2),'r','LineWidth',3); hold on
H2 = plot(tax(t),mean(bPower_AVHVH,2),'b','LineWidth',3); 
legend([H1,H2],{'AVL-VL','AVH-VH'},'Location','NorthWest')
xlim([0 0.5]); 
ylim([-0.8 0.8]); hline(0,':k')
xlabel('Time')
title('Beta Power','FontSize',14); ylabel('Power Vals','FontSize',14); 

saveas(gcf,'Figure_BetaTheta_RawVals_Conditions_OverTime.fig')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CORRELATE WEIGHT CHANGE WITH FREQ CHANGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

betaDiff = bPower_AVLVL - bPower_AVHVH;
betaDiff = betaDiff';
thetaDiff = tPower_AVLVL - tPower_AVHVH; 
thetaDiff = thetaDiff';

save('Beta','betaPower_AVLVL','betaPower_AVHVH','betaDiff'); 
save('Theta','thetaPower_AVLVL','thetaPower_AVHVH','thetaDiff');

% load in the perceptual weights 
load('W:/Project 1/log/EEG Data/allGroupBeh','awObs'); 
weightDiff = awObs(:,2)-awObs(:,1); 
for k = 1:7; 
    [r(k),p(k),stats{k}] = spearmanrank(weightDiff,betaDiff(:,k)); 
    [r2(k),p2(k),stats2{k}] = spearmanrank(weightDiff,thetaDiff(:,k));
end

% none 

%%
clear all; close all; clc 
load('Group_AVH.mat')
load('Group_AVL.mat'); 
load('Group_VH.mat'); 
load('Group_VL.mat'); 
load('Group_AUD.mat'); 

% Work out the diffs
AVL_VL = AVL-VL;
AVH_VH = AVH-VH; 

clearvars -except AVL_VL AVH_VH fax tax
% thetaE = [6 7 12 13 14 18 19 60]; 
thetaE = [12 13 14 18 19]; 
betaE = [7 8 9 12 13 14 18]; 
theta = find(fax>=4 & fax<=6); 
beta = find(fax>=10 & fax<=12); 
thetaT = find(tax>=0 & tax<=0.15);                        % shorter time window
betaT = find(tax>=0.3 & tax<=0.35); 


% POWER 
thetaLow = sq(AVL_VL(thetaE,theta,thetaT,:)); 
thetaHigh = sq(AVH_VH(thetaE,theta,thetaT,:)); 
betaLow = sq(AVL_VL(betaE,beta,betaT,:));
betaHigh = sq(AVH_VH(betaE,beta,betaT,:));

thetaDiff = [sq(mean(mean(mean(thetaHigh,1),2),3))  sq(mean(mean(mean(thetaLow,1),2),3))];
thetaDiff(:,3) = thetaDiff(:,2)-thetaDiff(:,1);    

% betaDiff = [sq(mean(betaHigh,1))'  sq(mean(betaLow,1))'];
% betaDiff(:,3) = betaDiff(:,2)-betaDiff(:,1);    

betaDiff = [sq(mean(mean(betaHigh,1),2))  sq(mean(mean(betaLow,1),2))];
betaDiff(:,3) = betaDiff(:,2)-betaDiff(:,1);  
% [AVH, AVL, AHL-AVH]

% BEHAVIOUR
load('W:/Project 1/log/EEG Data/groupData','BWeights','TAX'); 
thetaTB = find(TAX>=0 & TAX<=0.15); 
betaTB = find(TAX>=0.25 & TAX<=0.3); 
% 18:25 is sig

thetaWeights = sq(mean(BWeights(:,thetaTB,[1 3]),2));
betaWeights = sq(mean(BWeights(:,betaTB,[1 3]),2));

% % betaWeights = [BWeights(:,betaTB,1) BWeights(:,betaTB,3) BWeights(:,betaTB,2) BWeights(:,betaTB,4)];
% thetaWeights = [sq(mean(BWeights(:,thetaTB,1),2)) sq(mean(BWeights(:,thetaTB,3),2)) sq(mean(BWeights(:,thetaTB,2),2)) sq(mean(BWeights(:,thetaTB,4),2))]; 
% betaWeights = [sq(mean(BWeights(:,betaTB,1),2)) sq(mean(BWeights(:,betaTB,3),2)) sq(mean(BWeights(:,betaTB,2),2)) sq(mean(BWeights(:,betaTB,4),2))]; 


close
subplot 321; 
for s = 1:20; 
    plot(1:2,thetaWeights(s,1:2),'-o','Color',[0.7 0.7 0.7]); hold on 
%     plot(3:4,thetaWeights(s,3:4),'-o','Color',[0.7 0.7 0.7]); hold on
end
hold on; 
boxplot(thetaWeights(:,1:2),'Colors','k');
set(gca,'XTick',1:4,'Xticklabel',{'AVH','AVL'})
% set(gca,'XTick',1:4,'Xticklabel',{'AVH A','AVL A','AVH V','AVL V'})
title('behavioural regression weights, time 0:0.15, (theta)')
ylim([0 10])
ylabel('b weight')

% plot the behavour 
subplot 322; 
for s = 1:20; 
    plot(1:2,betaWeights(s,1:2),'-o','Color',[0.7 0.7 0.7]); hold on 
%     plot(3:4,betaWeights(s,3:4),'-o','Color',[0.7 0.7 0.7]); hold on
end
hold on; 
boxplot(betaWeights(:,1:2),'Colors','k');
set(gca,'XTick',1:2,'Xticklabel',{'AVH','AVL'})
% set(gca,'XTick',1:4,'Xticklabel',{'AVH A','AVL A','AVH V','AVL V'})
set(gca,'XTick',1:4,'Xticklabel',{'AVH A','AVL A','AVH V','AVL V'})
title('behavioural regression weights, time 0.3, (beta)')
ylim([0 35]); ylabel('b weight')

subplot 323; 
for s = 1:20; 
    plot(1:2,thetaDiff(s,1:2),'-o','Color',[0.7 0.7 0.7]); hold on 
end
hold on; 
boxplot(thetaDiff(:,1:2),'Colors','k');
set(gca,'XTick',1:4,'Xticklabel',{'AVH','AVL'})
title('theta power, time 0:0.15, (theta)')
ylabel('power'); ylim([-0.7 1.3]); hline(0); 


subplot 324
% plot the behavour 
for s = 1:20; 
    plot(1:2,betaDiff(s,1:2),'-o','Color',[0.7 0.7 0.7]); hold on 
end
hold on; 
boxplot(betaDiff(:,1:2),'Colors','k');
set(gca,'XTick',1:4,'Xticklabel',{'AVH','AVL'})
title('beta power, time 0:0.15, (beta)')
ylabel('power'); ylim([-0.7 1.3]); hline(0); 



% figure
% scatterplot
thetaD = thetaWeights(:,2)-thetaWeights(:,1); 
thetaPD = thetaDiff(:,3); 
betaD = betaWeights(:,2) - betaWeights(:,1); 
betaPD = betaDiff(:,3); 

% clf
subplot 325
scatterplot(thetaD,thetaPD);
axis([-2 4 -1 1 ])
% axis([-0.65 0.65 -0.65 0.65])
vline(0); hline(0); 
% vline(0,'Color',[0.7 0.7 0.7]); hline(0,'Color',[0.7 0.7 0.7])
xlabel('Behavioural Regression Weight Diff (AVL-AVH)'); 
ylabel('Theta Power Weight Diff (AVL-AVH)'); 
title('time window 0:150ms')

subplot 326
scatterplot(betaD,betaPD)
axis([-4.5 8.5 -1 1 ])
% vline(0,'Color',[0.7 0.7 0.7]); hline(0,'Color',[0.7 0.7 0.7])
vline(0); hline(0); 
xlabel('Behavioural Regression Weight Diff (AVL-AVH)'); 
ylabel('Beta Power Weight Diff (AVL-AVH)'); 
title('time window 300ms')
% saveas(gcf,'Figure3_NEW.fig')

% STATS

% stats_beh_theta = mes(thetaWeights(:,1),thetaWeights(:,2),'hedgesg','isDep',1);     % behavioural weights, 0:150ms
% stats_beh_theta_vis = mes(thetaWeights(:,3),thetaWeights(:,4),'hedgesg','isDep',1); 
% stats_power_theta = mes(thetaDiff(:,1),thetaDiff(:,2),'hedgesg','isDep',1);
% 
% stats_beh_beta = mes(betaWeights(:,1),betaWeights(:,2),'hedgesg','isDep',1);     % behavioural weights, 0:150ms
% stats_beh_beta_vis = mes(betaWeights(:,3),betaWeights(:,4),'hedgesg','isDep',1); 
% stats_power_beta = mes(betaDiff(:,1),betaDiff(:,2),'hedgesg','isDep',1);

%


bweightsDiff = betaWeights(:,2)-betaWeights(:,1); 
tweightsDiff = thetaWeights(:,2)-thetaWeights(:,1); 


[r,p,t] = spearmanrank(bweightsDiff,betaDiff(:,3))
[r,p,t] = spearmanrank(tweightsDiff,thetaDiff(:,3))





%%















% %%
% clear all; close all; clc 
% load('Group_AVH.mat')
% load('Group_AVL.mat'); 
% load('Group_VH.mat'); 
% load('Group_VL.mat'); 
% load('Group_AUD.mat'); 
% 
% % Work out the diffs
% AVL_VL = AVL-VL;
% AVH_VH = AVH-VH; 
% 
% clearvars -except AVL_VL AVH_VH fax tax
% e = 1:64; 
% theta = find(fax>=4 & fax<=6); 
% beta = find(fax>=10 & fax<=12); 
% thetaT = find(tax>=0 & tax<=0.15);                        % shorter time window
% betaT = find(tax==0.3); 
% 
% 
% % POWER 
% thetaLow = sq(AVL_VL(1:64,theta,thetaT,:)); 
% thetaHigh = sq(AVH_VH(1:64,theta,thetaT,:)); 
% betaLow = sq(AVL_VL(1:64,beta,betaT,:));
% betaHigh = sq(AVH_VH(1:64,beta,betaT,:));
% 
% thetaDiff = [sq(mean(mean(mean(thetaHigh,1),2),3))  sq(mean(mean(mean(thetaLow,1),2),3))];
% thetaDiff(:,3) = thetaDiff(:,2)-thetaDiff(:,1);    
% 
% betaDiff = [sq(mean(betaHigh,1))'  sq(mean(betaLow,1))'];
% betaDiff(:,3) = betaDiff(:,2)-betaDiff(:,1);    
% 
% % BEHAVIOUR
% load('W:/Project 1/log/EEG Data/groupData','BWeights','TAX'); 
% thetaTB = find(TAX>=0 & TAX<=0.15); 
% betaTB = find(TAX==0.3); 
% betaWeights = [BWeights(:,betaTB,1) BWeights(:,betaTB,3) BWeights(:,betaTB,2) BWeights(:,betaTB,4)];
% thetaWeights = [sq(mean(BWeights(:,thetaTB,1),2)) sq(mean(BWeights(:,thetaTB,3),2)) sq(mean(BWeights(:,thetaTB,2),2)) sq(mean(BWeights(:,thetaTB,4),2))]; 
% 
% %%
% clf
% subplot 221; 
% for s = 1:20; 
%     plot(1:2,thetaWeights(s,1:2),'-o','Color',[0.7 0.7 0.7]); hold on 
% %     plot(3:4,thetaWeights(s,3:4),'-o','Color',[0.7 0.7 0.7]); hold on
% end
% hold on; 
% boxplot(thetaWeights(:,1:2),'Colors','k');
% set(gca,'XTick',1:4,'Xticklabel',{'AVH','AVL'})
% % set(gca,'XTick',1:4,'Xticklabel',{'AVH A','AVL A','AVH V','AVL V'})
% title('behavioural regression weights, time 0:0.15, (theta)')
% 
% % plot the behavour 
% subplot 222; 
% for s = 1:20; 
%     plot(1:2,betaWeights(s,1:2),'-o','Color',[0.7 0.7 0.7]); hold on 
% %     plot(3:4,betaWeights(s,3:4),'-o','Color',[0.7 0.7 0.7]); hold on
% end
% hold on; 
% boxplot(betaWeights(:,1:2),'Colors','k');
% set(gca,'XTick',1:4,'Xticklabel',{'AVH','AVL'})
% % set(gca,'XTick',1:4,'Xticklabel',{'AVH A','AVL A','AVH V','AVL V'})
% set(gca,'XTick',1:4,'Xticklabel',{'AVH A','AVL A','AVH V','AVL V'})
% title('behavioural regression weights, time 0.3, (beta)')
% 
% subplot 223; 
% for s = 1:20; 
%     plot(1:2,thetaDiff(s,1:2),'-o','Color',[0.7 0.7 0.7]); hold on 
% end
% hold on; 
% boxplot(thetaDiff(:,1:2),'Colors','k');
% set(gca,'XTick',1:4,'Xticklabel',{'AVH','AVL'})
% title('theta power, time 0:0.15, (theta)')
% 
% subplot 224
% % plot the behavour 
% for s = 1:20; 
%     plot(1:2,betaDiff(s,1:2),'-o','Color',[0.7 0.7 0.7]); hold on 
% end
% hold on; 
% boxplot(betaDiff(:,1:2),'Colors','k');
% set(gca,'XTick',1:4,'Xticklabel',{'AVH','AVL'})
% title('beta power, time 0:0.15, (beta)')
% 
% %% scatterplot
% 
% thetaD = thetaWeights(:,2)-thetaWeights(:,1); 
% thetaPD = thetaDiff(:,3); 
% betaD = betaWeights(:,2) - betaWeights(:,1); 
% betaPD = betaDiff(:,3); 
% 
% clf
% subplot 121
% scatterplot(thetaD,thetaPD);
% axis([-2 4 -0.65 0.65 ])
% % axis([-0.65 0.65 -0.65 0.65])
% vline(0,'Color',[0.7 0.7 0.7]); hline(0,'Color',[0.7 0.7 0.7])
% xlabel('Behavioural Regression Weight Diff (AVL-AVH)'); 
% ylabel('Theta Power Weight Diff (AVL-AVH)'); 
% title('time window 0:150ms')
% 
% subplot 122;
% scatterplot(betaD,betaPD)
% axis([-4.5 8.5 -0.65 0.65 ])
% vline(0,'Color',[0.7 0.7 0.7]); hline(0,'Color',[0.7 0.7 0.7])
% xlabel('Behavioural Regression Weight Diff (AVL-AVH)'); 
% ylabel('Beta Power Weight Diff (AVL-AVH)'); 
% title('time window 300ms')
% 
% %% STATS
% 
% stats_beh_theta = mes(thetaWeights(:,1),thetaWeights(:,2),'hedgesg','isDep',1);     % behavioural weights, 0:150ms
% stats_beh_theta_vis = mes(thetaWeights(:,3),thetaWeights(:,4),'hedgesg','isDep',1); 
% stats_power_theta = mes(thetaDiff(:,1),thetaDiff(:,2),'hedgesg','isDep',1);
% 
% stats_beh_beta = mes(betaWeights(:,1),betaWeights(:,2),'hedgesg','isDep',1);     % behavioural weights, 0:150ms
% stats_beh_beta_vis = mes(betaWeights(:,3),betaWeights(:,4),'hedgesg','isDep',1); 
% stats_power_beta = mes(betaDiff(:,1),betaDiff(:,2),'hedgesg','isDep',1);
% 
% 
% 
% %%
% % 
% % %%
% % 
% % % differences in weights 
% % bweightDiff = [BWeights(:,betaTB,1) BWeights(:,betaTB,3)];
% % bweightDiff(:,3) = bweightDiff(:,2)-bweightDiff(:,1);
% % 
% % tweightDiff =  [sq(mean(BWeights(:,thetaTB,1),2)) sq(mean(BWeights(:,thetaTB,3),2))];
% % tweightDiff(:,3) = tweightDiff(:,2) - tweightDiff(:,1); 
% % 
% % % load in the perceptual weights 
% % load('W:/Project 1/log/EEG Data/groupPsyWeights','awObs'); 
% % awObs(:,3) = awObs(:,2)-awObs(:,1); 
% % 
% % [betaR,betaP] = corr(bweightDiff(:,3),betaDiff(:,3));
% % [thetaR,thetaP] = corr(tweightDiff(:,3),thetaDiff(:,3));
% % 
% % clf
% % subplot 121; 
% % scatterplot(awObs(:,3),betaDiff(:,3)); 
% % axis([-0.65 0.65 -0.65 0.65])
% % % line([-0.6,1],[-0.6,1],'Color',[0.7 0.7 0.7])
% % vline(0,'Color',[0.7 0.7 0.7]); hline(0,'Color',[0.7 0.7 0.7])
% % xlabel('WobsAVL - WobsAVH'); ylabel('Power AVL - Power AVH')
% % subplot 122; 
% % scatterplot(awObs(:,3),thetaDiff(:,3));
% % axis([-0.65 0.65 -0.65 0.65])
% % % line([-0.6,1],[-0.6,1],'Color',[0.7 0.7 0.7])
% % vline(0,'Color',[0.7 0.7 0.7]); hline(0,'Color',[0.7 0.7 0.7])
% % xlabel('WobsAVL - WobsAVH'); ylabel('Power AVL - Power AVH')
% % 
% % 
% % subplot 121
% % boxplot(
% % 
% % % 
% % 
% % 
% %     [r(k),p(k),stats{k}] = spearmanrank(weightDiff,betaDiff(:,k)); 
% %     [r2(k),p2(k),stats2{k}] = spearmanrank(weightDiff,thetaDiff(:,k));
% % 
% % for s =1:20; 
% %     plot(awObs(s,1:2),'-o','Color',[0.7 0.7 0.7]); hold on
% % end
% % boxplot(awObs(:,1:2))
% % subplot 311
% % boxplot(awObs(:,1:2))
% % subplot 312
% % boxplot(betaDiff(:,1:2))
% % subplot 313
% % boxplot(thetaDiff(:,1:2))
% % 
% % 
% % %% 
% % clf
% % subplot 121; 
% % scatterplot(awObs(:,3),betaDiff(:,3)); 
% % axis([-0.65 0.65 -0.65 0.65])
% % % line([-0.6,1],[-0.6,1],'Color',[0.7 0.7 0.7])
% % vline(0,'Color',[0.7 0.7 0.7]); hline(0,'Color',[0.7 0.7 0.7])
% % xlabel('WobsAVL - WobsAVH'); ylabel('Power AVL - Power AVH')
% % subplot 122; 
% % scatterplot(awObs(:,3),thetaDiff(:,3));
% % axis([-0.65 0.65 -0.65 0.65])
% % % line([-0.6,1],[-0.6,1],'Color',[0.7 0.7 0.7])
% % vline(0,'Color',[0.7 0.7 0.7]); hline(0,'Color',[0.7 0.7 0.7])
% % xlabel('WobsAVL - WobsAVH'); ylabel('Power AVL - Power AVH')
