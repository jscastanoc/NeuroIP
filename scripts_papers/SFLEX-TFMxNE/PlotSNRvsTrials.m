% Plot SNR vs Trials

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SNR vs Trials %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 19 Aug 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization and creation of the structures used in some functions
close all; clc; clear

addpath('../../external/source_toolbox/haufe/')
addpath('../../external/source_toolbox/nolte/')
addpath('../../external/source_toolbox/simulations/')
rng('default')
rng('shuffle');

warning off

load clab_example;
load clab_10_10;
clab = clab_10_10;
% data_name = 'icbm152b_sym';
data_name = 'montreal';

sa = prepare_sourceanalysis(clab, data_name);

temp = sa.V_cortex_coarse;
L = nip_translf(temp); % Leadfield matrix
L = L(find(ismember(clab_example,clab)),:);
clear temp

% cfg.cortex = sa.cortex10K;
cfg.L = L;
cfg.fs = 200; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:2; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa;


%% Generation of the Morlet wavelet %%


clear source;

%% Simulation of the EEG (Ntrial Averages)
Ntrials = [5 20 50 100 250];
snr_meas = 0;
snr_bio = [-15 -5 5];
Nspurious = 500;
snr = [];
t = model.t;
L = model.L;
% cortex = model.cortex;
clear model;
save_trials = true;
n = 1;
lspec = {'-k','--k','-.k'};


snr_fig = figure('Units','normalized','position',[0.1 0.1 0.2 0.3]);
for i = snr_bio
    for j = 1:50
        snr = [];
        for k = Ntrials;
            dir = 'D:/Datasets/sim_trials_final/1/';
            file_name = strcat(dir,'Exp',num2str(j),'Ntrials',...
                num2str(k),'BioNoise',num2str(i),'.mat');
            load(file_name);
            y = mgjob.results{1};
            Jclean = mgjob.results{2};
            actidx = mgjob.results{3};
            y_clean = L*Jclean;
            noise = y-y_clean;
            snr = [snr, 20*log10(norm(y_clean)/norm(noise))];
        end
        snr_vec(j,:) = snr;
    end
    figure(snr_fig)
    inter_vec = (5:5:250);
    plot(inter_vec,interp1(Ntrials,mean(snr_vec),inter_vec,'cubic'),lspec{n})
    hold on
    n =n+1;
end
pause
% save('snr_',num2str(k),'trials.mat')

