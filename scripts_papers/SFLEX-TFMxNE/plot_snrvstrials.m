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
data_name = 'icbm152b_sym';
% data_name = 'montreal';

sa = prepare_sourceanalysis(clab, data_name);

temp = sa.V_cortex10K;
L = nip_translf(temp); % Leadfield matrix
L = L(find(ismember(clab_example,clab)),:);
clear temp

cfg.cortex = sa.cortex10K;
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
Ntrials = [5 50 100 250 500];
snr_meas = 0;
snr_bio = -5;
Nspurious = 500;
snr = [];

for j = 1:20
    j
    for i = Ntrials
    [model.y, Jclean] = nip_simtrials(model.L, model.cortex.vc, ...
            source_act, model.t, Nspurious , i, snr_meas, snr_bio);
        y_clean = model.L*Jclean;
        noise = model.y-y_clean;
        snr = [snr 20*log10(norm(y_clean)/norm(noise))];
    end
    snr_array(j,:) = snr;
    snr = [];
end
plot(Ntrials,mean(snr_array))
figure('Units','normalized','position',[0.1 0.1 0.2 0.3]);
plot(5:15:500,interp1(Ntrials,mean(snr_array),5:15:500,'cubic'),'k');
ylabel('SNR (dB)'); xlabel('# Trials'); ylim([-6 10])

    