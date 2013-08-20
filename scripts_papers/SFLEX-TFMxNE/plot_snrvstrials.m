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


clear source;

%% Simulation of the EEG (Ntrial Averages)
Ntrials = [5 50 100 250 500];
snr_meas = 0;
snr_bio = [-15 -5 5];
Nspurious = 500;
snr = [];
t = model.t;
L = model.L;
cortex = model.cortex;
clear model;
save_trials = true;
for k = 1:length(snr_bio)
    for j = 1:50       
        phase_shift = [0.5 1.5] + 0.15*randn(1,2); % Phase shift for the sources (in seconds)
        Nact = length(phase_shift);
        fc_wl = [5 5] + 1*randn(1,2); % Central frequency for the wavelet
%         for i = 1:Nact
%             f0 = fc_wl(i);
%             
%             % Normalization terms for the wavelet
%             sigma_f = f0/7;
%             sigma_t = 1/(2*pi*sigma_f);
%             
%             % "source" contains the time courses of active sources
%             source(i,:) =  real(exp(2*1i*pi*f0*model.t).*...
%                 exp((-(model.t-phase_shift(i)).^2)/(2*sigma_t^2)));
%         end
%         source_act = source;
        source_act = [aux_wavelet(fc_wl(1),t,phase_shift(1)); aux_wavelet(fc_wl(2),t,phase_shift(2))];
        snr = [];
        for i = Ntrials
            [y, Jclean, actidx] = nip_simtrials(L, cortex.vc, ...
                source_act, t, Nspurious , i, snr_meas, snr_bio(k));
            y_clean = L*Jclean;
            noise = y-y_clean;
%             y = model.y;
            snr = [snr 20*log10(norm(y_clean)/norm(noise))];
            if save_trials     
                dir = 'D:/Datasets/sim_trials/';
                Jclean = sparse(Jclean);
                file_name = strcat(dir,'Exp',num2str(j),'Ntrials',...
                    num2str(i),'BioNoise',num2str(snr_bio(k)),'.mat');
                save(file_name,'y','y_clean','Jclean','actidx');
            end
        end
        snr_array{k}(j,:) = snr;
    end    
end

save('snrvstrials.mat','snr_array')

plot(Ntrials,mean(snr_array{2}))
figure('Units','normalized','position',[0.1 0.1 0.2 0.3]);
plot(5:15:500,interp1(Ntrials,mean(snr_array),5:15:500,'cubic'),'k');
ylabel('SNR (dB)'); xlabel('# Trials'); ylim([-6 10])

    