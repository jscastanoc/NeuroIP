% Script BASIC
clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add data/ and external libraries to the path%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% External libraries needed for STOUT:
% LTFAT
% nip_init();

% ltfatstart; % Only need to be done once per matlab session


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated Brain Activity %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Forward problem data.
load(strcat('data/sa_montreal.mat'))

cfg.L = nip_translf(sa.V_cortex_coarse);
cfg.cortex = sa.cortex_coarse;
cfg.fs = 120; % Sample frequency
cfg.t = 0:1/cfg.fs:1; % Time vector

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

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


% Simulate 2 active dipoles
% [J, ~] = nip_simulate_activity(model.cortex.vc,[30 -20 30], act, randn(1,3), model.t);
[J, ~] = nip_simulate_activity(model.cortex.vc, [30 -20 30;-30 20 30] , act,  ones(size(act,1),3), model.t, struct('sample_all',0));


% Simulate "smooth activity" (simulation is spatial low pass filtered
fuzzy = nip_fuzzy_sources(model.cortex,1.5);
index = (1:3:model.Nd);

for i = 0:2
    J(index+i,:) = fuzzy*J(index+i,:);
end
% J contains the final simulated activity


% Obtain pseudo EEG corresponding to the simulation
clean_y = model.L*J;

% Add noise at sensor level
snr = 10;
model.y = nip_addnoise(clean_y, snr);

transM = (1+(1/model.Nc))*eye(model.Nc)-(1/model.Nc)*ones(model.Nc,model.Nc);

model.y = transM*model.y;
model.L = transM*model.L;

% Depth compensation
depth = 'Lnorm'; % it can be none, Lnorm or sLORETA-based depth compensation
switch depth
    case 'none'
        L = model.L;
        Winv = [];
    case 'Lnorm'
        gamma = 0.7; % How strong the depth compensation is?
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;
end
clear extras;

%% SOLUTION %%


sigma = 1; % Width of the gaussian bells or cortical blobs used as spatial dictionary

% This function can be used to create spatial dictionaries using a forward
% model taking into account only the cortex surface (model.cortex has the
% 3d graph representing the cortical surface)
B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));

% Normalize the spatial basis functions
B = nip_blobnorm(B,'norm',2);



% Options for the inversion
% The ratio between spatial_reg and temp_reg depends on the snr of the EEG
% However, in general a ratio of 1:3 should work ok.
spatial_reg = 5 ; % Sparsity in the spatial domain
temp_reg =  1; % Sparsity in the time-frequency domain

% Set regularization parameters to get an ideal goodness of fit
gof = norm(model.y - model.L*J, 'fro')/norm(model.y, 'fro')

a = 10;  %  Time shift for the Short Time Fourier Transform (STFT).
m = 100; %Frequency bins for the STFT.
lipschitz = [];
% By setting 'optimgof' to true, the regularization parameters will be
% modified to get the desired gof. If false, then the user-select reg.
% parameters are used for the solution
[J_est, extras] = nip_stout(model.y, L, B,'optimgof',true,...
    'sreg',spatial_reg,'treg',temp_reg,'gof', gof, 'a',a ,'m',m,...
    'lipschitz', lipschitz,'Winv',Winv);

resnorm = norm(model.y-model.L*J_est, 'fro')/norm(model.y, 'fro')

%% Visualization %%
%%% Simulation
%Temporal
figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
plot(model.t,J')
xlabel('Time')
ylabel('Amplitude')

% Spatial
figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
nip_reconstruction3d(model.cortex, sqrt(sum(J.^2,2)), struct('axes',gca)); 



%%% Reconstruction %%%
%Temporal
figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
plot(model.t,J_est')
xlabel('Time')
ylabel('Amplitude')

% Spatial
figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
nip_reconstruction3d(model.cortex, sqrt(sum(J_est.^2,2)),  struct('axes',gca)); 

