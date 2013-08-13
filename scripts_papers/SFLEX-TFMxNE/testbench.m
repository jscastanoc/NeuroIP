%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Testbench for LORETA SFLEX TFMxNE and Proposed %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 13 Aug 2013           %
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
% sa = prepare_sourceanalysis(clab, 'icbm152b_sym');
sa = prepare_sourceanalysis(clab, 'montreal');

temp = sa.V_cortex;
L = reshape(permute(temp,[1 3 2]), size(temp,1), size(temp,2)*3); % Leadfield matrix
clear temp

cfg.cortex = sa.cortex;
cfg.L = L;
cfg.fs = 200; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:2; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa;


%% Generation of the Morlet wavelet %%

phase_shift = [0.5 1.5] ; % Phase shift for the sources (in seconds)
Nact = length(phase_shift);
fc_wl = 5; % Central frequency for the wavelet
for i = 1:Nact
    f0 = fc_wl;
    
    % Normalization terms for the wavelet
    sigma_f = f0/7;
    sigma_t = 1/(2*pi*sigma_f);
    
    % "source" contains the time courses of active sources
    source(i,:) =  real(exp(2*1i*pi*f0*model.t).*...
        exp((-(model.t-phase_shift(i)).^2)/(2*sigma_t^2)));
end
source_act = source;
clear source;

%% Simulation of the EEG (Ntrial Averages)
Ntrials = 50;
snr_meas = 10;
snr_bio = 0;
Nspurious = 200;
[y, Jclean] = nip_simtrials(model.L, model.cortex, ...
        source_act, model.t, Nspurious , Ntrials, snr_meas, snr_bio);
    
figure;
nip_reconstruction3d(model.cortex,sqrt(sum(Jclean.^2,2)),gca);
    
fuzzy = nip_fuzzy_sources(model.cortex,0.1);
    