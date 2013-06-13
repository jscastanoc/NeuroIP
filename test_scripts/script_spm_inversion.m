% MNE with ReML (free energy as cost function)

% Script BASIC
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
cfg.fs = 500; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:0.5; % Vector de tiempo

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec


% Extraer el laplaciano espacial (ahí se codifica la información de vecinos
[Laplacian, ~] = nip_neighbor_mat(model.cortex);

act = sin(2*pi*10*model.t); % Actividad a simular
% Simular actividad "act" en el dipolo más cercano a las coordenadas [30
% -20 30]
J = nip_simulate_activity(model.cortex,Laplacian, [30 -20 30], ...
        act,model.t);

% Hay dispersion de la actividad, entre más grande el número, más dispersa es la actividad    
fuzzy = nip_fuzzy_sources(model.cortex,1);
J = fuzzy*J; % J simulado FINAL

% Obtener el eeg correspondiente a la simulacion
clean_y = model.L*J;

% Anadir ruido
snr = 10;
model.y = nip_addnoise(clean_y, snr);



%%%%%%%%%%%%%%
% Estimacion %
%%%%%%%%%%%%%%
% Estimar actividad (en este caso LORETA por que se usa el laplaciano
% espacial para hallar la matriz de covarianza
Q = inv(Laplacian'*Laplacian); %Matriz de covarianza apriori
Qe = eye(model.Nc);
[J_est, extras_spm] = nip_spm_inversion(model.y, model.L, Q,Qe);
Q = diag(nip_lcmv(model.y,model.L));
[J_est2, extras_spm2] = nip_spm_inversion(model.y, model.L, Q,Qe);
% [J_lor, extras_lor] = nip_loreta(model.y, model.L, Q);
% J_est = nip_sloreta(model.y,model.L);

figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
nip_reconstruction3d(model.cortex,sqrt(sum(J_est.^2,2)),gca);
title('ReML')

figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
nip_reconstruction3d(model.cortex,sqrt(sum(J_est2.^2,2)),gca);
title('GCV')