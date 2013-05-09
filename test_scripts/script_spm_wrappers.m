% Script to test ARD, GS and MSP

clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 4000; % 1000, 2000, 4000, 8000

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('data/montreal',num2str(Nd),'_full.mat'))

cfg.L = L;
cfg.cortex = cortex_mesh;
cfg.fs = 200; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:1; % Vector de tiempo

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec


% Extraer el laplaciano espacial (ahí se codifica la información de vecinos
[Laplacian, ~] = nip_neighbor_mat(model.cortex);

act(1,:) = sin(2*pi*10*model.t); % Actividad a simular
act(2,:) = sin(2*pi*10*model.t);
J = nip_simulate_activity(model.cortex,Laplacian, [30 -20 30;-30 -20 30], ...
        act,model.t);

% Hay dispersion de la actividad, entre más grande el número, más dispersa es la actividad    
fuzzy = nip_fuzzy_sources(model.cortex,0.5);
J = fuzzy*J; % J simulado FINAL

figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
subplot(1,3,1)
nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),gca);


% Obtener el eeg correspondiente a la simulacion
clean_y = model.L*J;

% Añadir ruido
snr = -10;
model.y = nip_addnoise(clean_y, snr);


% Generacion de parches
Np = 256; % Numero de parches
fuzzy = nip_fuzzy_sources(model.cortex,1.5);
fuzzy = fuzzy(:,round(linspace(1,model.Nd,Np)));

% Priors con GS
h = nip_spm_solvers(model.y,model.L,fuzzy,eye(model.Nc),'GS');
subplot(1,3,2)
nip_reconstruction3d(model.cortex,fuzzy*h,gca);

% Priors con ARD
h = nip_spm_solvers(model.y,model.L,fuzzy,eye(model.Nc),'ARD');
subplot(1,3,3)
nip_reconstruction3d(model.cortex,fuzzy*h,gca);