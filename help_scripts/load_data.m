% Script to load the data I'm typically using

% Load lead field, cortex mesh, labels etc...
load(strcat('data/sa_montreal.mat'))
load clab_full % Labels for 118 channels
clab_full = labels;
load('clab_10_10_34c.mat'); % Labels for 59 channels under 10 10 protocol
% clab = clab_10_10;
clab = labels;
% data_name = 'montreal';
% sa = prepare_sourceanalysis(clab, data_name);

% Select rows of the lead field according of the EEG standard to be used
cfg.L = nip_translf(sa.V_cortex_coarse);
cfg.L = cfg.L(find(ismember(clab_full,clab)),:);
cfg.cortex = sa.cortex_coarse;
cfg.clab = clab;

% Set time information
cfg.fs = 120; % Sample rate (in Hz)
cfg.t = 0:1/cfg.fs:1.5; % Time vector

% model structure with Nd, Nc and Nt according to what we defined above
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec