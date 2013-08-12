% Temporal Proyection with dct

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
act(2,:)  = exp(-(model.t-0.75).^2/0.05).*sin(2*pi*10*model.t);


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

figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),gca)

figure('Units','normalized','position',[0.2 0.2 0.14 0.3]);
nip_sens_meas(model.y,model.labels,gca)
title('Raw data')

y_dct = dct(model.y')';
figure('Units','normalized','position',[0.2 0.2 0.14 0.3]);
nip_sens_meas(y_dct,model.labels,gca)


yty = y_dct'*y_dct;
figure('Units','normalized','position',[0.2 0.2 0.3 0.3]);
imagesc(yty)
axis image

[Ur, S] = svd(yty);
Nr = 2;
temp_proj = Ur(:,1:Nr);

figure('Units','normalized','position',[0.2 0.2 0.3 0.3]);
plot(temp_proj)
title('Temporal Projector')

yp = model.y*idct(temp_proj);

figure('Units','normalized','position',[0.2 0.2 0.3 0.3]);
plot(yp')
title('Projected Data')


figure('Units','normalized','position',[0.2 0.2 0.3 0.3]);
Q = diag(nip_lcmv(yp,model.L));
[J_est,~] = nip_loreta(yp,model.L,Q);
nip_reconstruction3d(model.cortex,sqrt(sum(J_est.^2,2)),gca);
title('Projected Results')

figure('Units','normalized','position',[0.2 0.2 0.3 0.3]);
Q = diag(nip_lcmv(model.y,model.L));
[J_est,~] = nip_loreta(model.y,model.L,Q);
nip_reconstruction3d(model.cortex,sqrt(sum(J_est.^2,2)),gca);
title('Raw Results')

[U S V] = svd(model.y'*model.y);
yp_juana = model.y*U(:,1:Nr);
figure('Units','normalized','position',[0.2 0.2 0.3 0.3]);
Q = diag(nip_lcmv(yp_juana,model.L));
[J_est,~] = nip_loreta(yp_juana,model.L,Q);
nip_reconstruction3d(model.cortex,sqrt(sum(J_est.^2,2)),gca);
title('Juana Results')




