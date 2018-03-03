%% Example to analyze *.nirs data
%
%Last Update: 2018-03-03
%
%By: Guilherme A. Zimeo Morais
%Contact: gazmorais@gmail.com

%% Set trigger markers in the *.nirs files
    %Note: please make sure it is working properly
    rest_dur = 20; %duration of rest blocks
    task_dur = 20; %duration of task blocks
    pathname = cd; %current directory as ref
    overwrt = 1; %flag to overwrite original files
    InputTriggers(rest_dur, task_dur, pathname, overwrt);

%% Load File(s)
    %Please organize files according to the hierarchy: \group\subject\.nirs
    %e.g. ...\G1\S1\Abc.nirs , ...\G1\S2\Xyz.nirs, ...\G2\S1\Pqr.nirs , etc
    %Note: you may use other hierarchies, but adapt 'folderHierarchy' below
    rootFolder = cd;
    folderHierarchy = {'group', 'subject'}; %hierarchy of the folders
    loadFunc = {@nirs.io.loadDotNirs};
    fileExt = {'.nirs'};
    raw = nirs.io.loadDirectory(rootFolder, folderHierarchy, loadFunc, fileExt);
    
%% Create Demographics Table
    demographics = nirs.createDemographicsTable(raw);

%% Set task duration (this follows the example provided by nirs-toolbox)
    task_dur = 20;
    raw_edit = nirs.design.change_stimulus_duration(raw,'stim_channel1',task_dur);
    
%% Remove "auxiliary" trigger marker
    job = nirs.modules.DiscardStims();
    job.listOfStims = {
                       'aux_channel1'; 
                      };
    raw_edit = job.run(raw_edit);
    
%% Truncate Time Series
    job = nirs.modules.TrimBaseline(); %truncate time series
    job.preBaseline  = 30; % prior to first onset
    job.postBaseline = 30; % after last onset + duration
    raw_edit = job.run(raw_edit);

%% Compute [Unfiltered] Hemodynamic Changes
    %for AR-IRLS, it is recommended *not* to apply any filter
    job = nirs.modules.OpticalDensity();
    job = nirs.modules.BeerLambertLaw(job);
    hb = job.run(raw_edit);

%% Individual Stats: AR-IRLS
    job = nirs.modules.AR_IRLS();
    job.verbose = true;
    job.basis = Dictionary();
    canonical = nirs.design.basis.Canonical();
    canonical.peakTime = 6; %as in original hrf
    job.basis('default') = canonical;
    SubjStats=job.run(hb);

%% Visualize Individual Results
    for i=1:length(SubjStats)
        pcor = 0.05/(length(SubjStats(i).probe.distances)/2);
        critV = ['p < ' num2str(pcor)]; %Bonferroni
        thres(i) = SubjStats(i).getCritT(critV); %get threshold t-value
        SubjStats(i).variables((abs(SubjStats(i).tstat)>thres(i)),:) %list channels
        SubjStats(i).draw('tstat', [-10 10], critV); %draw 2D results
    end

%% Group Stats
    job = nirs.modules.MixedEffects();
    job.formula = 'beta ~ -1 + group:cond + (1|subject)';
    job.dummyCoding = 'full';
    GroupStats = job.run(SubjStats);

%% Visualize Group Results
    pcor = 0.05/(length(GroupStats.probe.distances)/2);
    critV = ['p < ' num2str(pcor)]; %Bonferroni
    thres = GroupStats.getCritT(critV);
    disp(GroupStats.table()); %summary table of results
    GroupStats.variables((abs(GroupStats.tstat)>thres),:) %list channels 
    GroupStats.draw('tstat', [-10 10], critV); %draw 2D results

%% Apply band-pass filter for data visualization
    job = eeg.modules.BandPassFilter();
    job.do_downsample = 0;
    job.highpass = 0.01;
    job.lowpass = 0.2;
    job.keepdc = 1;
    hb_filt = job.run(hb);

%% Plot Block Average
    blen = 5; %baseline prior onset to correct for
    offset = 15; %offset after task duration to plot
    chnlst = []; %list of channels to plot (empty = all)
    keep = [1 1 0 1 1 0 0 1 0 1]; %channels to be included
    keep = repmat(keep, length(hb_filt), 1); %all subjects
    plotflag = [1 0 0]; %flags for plots (details in code)
    hbnum = [1 2]; %plot both HbO and HbR curves
    BlockAvg(hb_filt, blen, offset, chnlst, keep, plotflag, hbnum);