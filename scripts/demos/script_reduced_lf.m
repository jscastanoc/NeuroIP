% Script BASIC
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add data/ and external libraries to the path%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% External libraries needed for STOUT:
% LTFAT
nip_init();

% ltfatstart; % Only need to be done once per matlab session


%%%%%%%%%%%%%%%%%%%%%%%%%%%%Hannah
% Simulated Brain Activity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Forward problem data.
% load(strcat('data/sa_montreal.mat'))
% 
% cfg.L = nip_translf(sa.V_coarse);
% cfg.cortex.vc = sa.grid_coarse;
% cfg.fs = 120; % Sample frequency
% cfg.t = 0:1/cfg.fs:1; % Time vector
% 
% % Crea una estructura (model) con los datos cargados y declarados arriba
% model = nip_create_model(cfg);
% clear cfg L cortex_mesh eeg_std head elec
load_data_grid;

%% Generate time series for the active dipoles

 % Phase shift for the sources (in seconds)
phase_shift = [0.3 0.6];

% Get number of active dipoles
Nact = length(phase_shift);

% Central frequency for the wavelets
fc_wl = [7 12]; 

% Generate time series
for i = 1:Nact
    f0 = fc_wl(i);
    
    % Normalization terms for the wavelet
    sigma_f = f0/7;
    sigma_t = 1/(2*pi*sigma_f);
    
    % "source" contains the time courses of active sources
    act(i,:) =  real(exp(2*1i*pi*f0*model.t).*...
        exp((-(model.t-phase_shift(i)).^2)/(2*sigma_t^2)));
end

rng(1);
% Simulate 2 active dipoles
dir_sim = randn(2,3);
% [J, actidx] = nip_simulate_activity(model.cortex.vc,[-30 20 30], act, dir_sim, model.t,struct('sample_all',1));
[J, actidx] = nip_simulate_activity(model.cortex.vc, [30 -20 30;-30 20 30] , act,  dir_sim, model.t, struct('sample_all',1));


% Simulate "smooth activity" (simulation is spatial low pass filtered
fuzzy = nip_fuzzy_sources(model.cortex,0.015);
index = (1:3:model.Nd);

for i = 0:2
    J(index+i,:) = fuzzy*J(index+i,:);
end
% J contains the final simulated activity


% Obtain pseudo EEG corresponding to the simulation
clean_y = model.L*J;

% Add noise at sensor level
snr = 5;
model.y = nip_addnoise_bio(model.L,J,model.t,snr);

% transM = eye(model.Nc)-(1/model.Nc)*ones(model.Nc);
% model.y = transM*model.y;
% model.L = transM*model.L;

% Depth compensation
depth = 'Lnorm'; % it can be none, Lnorm or sLORETA-based depth compensation
switch depth
    case 'none'
        L = model.L;
        Winv = [];
    case 'Lnorm'
        gamma = 0.6;
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;
end
clear extras;

%% SOLUTION %%


sigma = 1.2; % Width of the gaussian bells or cortical blobs used as spatial dictionary

% This function can be used to create spatial dictionaries using a forward
% model taking into account only the cortex surface (model.cortex has the
% 3d graph representing the cortical surface)
B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));

% Normalize the spatial basis functions
B = nip_blobnorm(B,'norm',2);



% Options for the inversion
% The ratio between spatial_reg and temp_reg depends on the snr of the EEG
% However, in general a ratio of 1:3 should work ok.
spatial_reg =90 ; % Sparsity in the spatial domain
temp_reg = 1; % Sparsity in the time-frequency domain

% Set regularization parameters to get an ideal goodness of fit
resnorm = norm(model.y - model.L*J, 'fro')/norm(model.y, 'fro')

a = 8;  %  Time shift for the Short Time Fourier Transform (STFT).
m = 64; %Frequency bins for the STFT.
lipschitz = [];
% By setting 'optimres' to true, the regularization parameters will be
% modified to get the desired resnorm. If false, then the user-select reg.

% temp_reg = 0;

% load('/mnt/data/Master_Results/Datasets/simulated/montreal_grid_orig/3/Exp11Ntrials100BioNoise-5.mat')
% load('/mnt/data/Master_Results/Datasets/simulated/montreal_grid_orig/3/Exp6Ntrials100BioNoise-5.mat')
% J = full(Jclean);
% resnorm = gof;
% model.y = y;

[J_eststout, extras] = nip_stout_python(model.y, L, B,'optimres',true,...
    'sreg',spatial_reg,'treg',temp_reg,'resnorm', resnorm, 'tstep',a ,'wsize',m,...
    'lipschitz', lipschitz,'Winv',Winv);
er{1}= nip_all_errors(model.y,model.L,J_eststout,J,model.cortex,actidx);


[J_estsflex, extras] = nip_sflex(model.y, L, B, 'optimres',true,'regpar',10,'resnorm', resnorm,'Winv',Winv);
er{2} = nip_all_errors(model.y,model.L,J_estsflex,J,model.cortex,actidx);

[J_estsflexFISTA, extrasF] = nip_sflex_python(model.y, L,B , 'optimres',true, 'resnorm',resnorm,...
    'sreg',spatial_reg,'Winv',Winv);
% save('J_estsflexFISTA.mat','J_estsflexFISTA');
er{3}= nip_all_errors(model.y,model.L,J_estsflexFISTA,J,model.cortex,actidx);

% break;
% 
% [J_esttfmxne, extras] = nip_tfmxne_python(model.y, L, 'optimres',true, 'resnorm',resnorm,...
%     'sreg',spatial_reg,'treg',temp_reg, 'tstep',a ,'wsize',m,...
%     'lipschitz', lipschitz,'Winv',Winv);
% er{3}= nip_all_errors(model.y,model.L,J_esttfmxne,J,model.cortex,actidx);


% J_est = J_estsflex;
% J_est = J_esttfmxne;
% J_est = J_eststout;

% resnorm = norm(model.y-model.L*J_est, 'fro')/norm(model.y, 'fro');




%% Visualization %%

source_pars = struct('orientation', 'coronal', 'ncol', [3 4], 'dslice_shown', 2.5, 'mydipmarkersize', 1.5, 'mydiplinewidth', 2.5, 'trcut', 0.25);

figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
h = showmri_transp(sa.mri, source_pars, [sa.grid_coarse 100*sum(nip_energy(J),2)], [sa.grid_coarse(actidx, :) dir_sim]);
zlab = get(h.cb, 'ylabel');
set(zlab, 'String', 'k.A.', 'fontsize', 18) 
set(h.cb, 'fontsize', 18)

figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
h = showmri_transp(sa.mri, source_pars, [sa.grid_coarse 100*sum(nip_energy(J_estsflex),2)], [sa.grid_coarse(actidx, :) dir_sim]);
zlab = get(h.cb, 'ylabel');
set(zlab, 'String', 'k.A.', 'fontsize', 18) 
set(h.cb, 'fontsize', 18)

figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
h = showmri_transp(sa.mri, source_pars, [sa.grid_coarse 100*sum(nip_energy(J_estsflexFISTA),2)], [sa.grid_coarse(actidx, :) dir_sim]);
zlab = get(h.cb, 'ylabel');
set(zlab, 'String', 'k.A.', 'fontsize', 18) 
set(h.cb, 'fontsize', 18)

figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
h = showmri_transp(sa.mri, source_pars, [sa.grid_coarse 100*sum(nip_energy(J_eststout),2)], [sa.grid_coarse(actidx, :) dir_sim]);
zlab = get(h.cb, 'ylabel');
set(zlab, 'String', 'k.A.', 'fontsize', 18) 
set(h.cb, 'fontsize', 18)
%%% Simulation
%Temporal
% figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
% plot(model.t,J')
% xlabel('Time')
% ylabel('Amplitude')
% 
% % Spatial
% figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
% nip_reconstruction3d(model.cortex, sqrt(sum(J.^2,2)), struct('axes',gca)); 
% 
% 
% 
% %%% Reconstruction %%%
% %Temporal
% figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
% plot(model.t,J_est')
% xlabel('Time')
% ylabel('Amplitude')
% 
% % Spatial
% figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
% nip_reconstruction3d(model.cortex, sqrt(sum(J_est.^2,2)),  struct('axes',gca)); 