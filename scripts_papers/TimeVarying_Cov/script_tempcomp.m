% Script test nip_tempcomp

clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 2000; % 1000, 2000, 4000, 8000

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('data/montreal',num2str(Nd),'_10-10.mat'))

cfg.L = L;
cfg.cortex = cortex_mesh;
cfg.fs = 500; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:1; % Vector de tiempo

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
model.labels = elec.label;
clear cfg L cortex_mesh eeg_std head elec


% Extraer el laplaciano espacial (ahí se codifica la información de vecinos
[Laplacian, ~] = nip_neighbor_mat(model.cortex);


act(1,:) = exp(-(model.t-0.25).^2/0.05).*sin(2*pi*10*model.t);
act(2,:)  = exp(-(model.t-0.75).^2/0.05).*sin(2*pi*20*model.t);


% Simular actividad "act" en el dipolo más cercano a las coordenadas [30
% -20 30]
J = nip_simulate_activity(model.cortex,Laplacian, [30 -20 30; -30 -20 30], ...
        act,model.t);

% Hay dispersion de la actividad, entre más grande el número, más dispersa es la actividad    
fuzzy = nip_fuzzy_sources(model.cortex,1);
J = fuzzy*J; % J simulado FINAL

% Obtener el eeg correspondiente a la simulacion
clean_y = model.L*J;

% Añadir ruido
snr = -10;
model.y = nip_addnoise(clean_y, snr);


[y_proj, y_rec, Ur, Er] = nip_tempcomp(model.y, model.t, [1 50], 0.7)


figure('Units','normalized','position',[0.2 0.2 0.15 0.5])
subplot(3,1,1)
plot(model.y')
subplot(3,1,2)
plot(y_rec')
subplot(3,1,3)
plot(Ur)