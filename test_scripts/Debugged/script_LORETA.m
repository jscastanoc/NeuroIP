% Script BASIC
clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 2000; % 2000, 4000

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('../../data/nocons_or/montreal',num2str(Nd),'_full1shell.mat'))

cfg.L = lf;
cfg.cortex = cortex_mesh;
cfg.fs = 200; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:1; % Vector de tiempo

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec


% Extraer el laplaciano espacial (ahí se codifica la información de vecinos
[Laplacian] = nip_neighbor_mat(model.cortex);


act = sin(2*pi*10*model.t); % Actividad a simular

% Simular actividad "act" en el dipolo más cercano a las coordenadas [30
% -20 30]
[J, ~] = nip_simulate_activity(model.cortex.vertices,[30 -20 30], act, randn(1,3), model.t);

% Hay dispersion de la actividad, entre mas grande el numero, mas dispersa es la actividad    
fuzzy = nip_fuzzy_sources(model.cortex,0.1);
index = (1:3:model.Nd);

for i = 0:2
    J(index+i,:) = fuzzy*J(index+i,:); % J simulado FINAL
end

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
Q = eye(model.Nd); %Matriz de covarianza apriori
[J_est, extras] = nip_loreta(model.y, model.L, Q);
% J_est = nip_sloreta(model.y,model.L);


%%%%%%%%%%%%%%%%%
% Visualizacion %
%%%%%%%%%%%%%%%%%
%%% Visualizacion de la simulacion %%%
%Temporal
figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
plot(model.t,J')
xlabel('Time')
ylabel('Amplitude')

% Espacial en un instante de tiempo dado
% t_0 = 25ms
t_0 = 0.025*model.fs;
figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
nip_reconstruction3d(model.cortex, J(:,t_0), gca);

%%% Visualizacion de la reconstruccion %%%
%Temporal
figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
plot(model.t,J_est')
xlabel('Time')
ylabel('Amplitude')

% Espacial en un instante de tiempo dado
% t_0 = 25ms
t_0 = 0.025*model.fs;
figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
nip_reconstruction3d(model.cortex, J_est(:,t_0), gca);
