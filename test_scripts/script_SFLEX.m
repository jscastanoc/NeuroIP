% Script SFLEX

clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 4000; % 1000, 2000, 4000, 8000

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('data/montreal',num2str(Nd),'_10-10.mat'))

cfg.L = L;
cfg.cortex = cortex_mesh;
cfg.fs = 500; % Frecuencia de muestreo para la simulacion
cfg.t = 0; % Vector de tiempo
cfg.elec = elec;

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

% act(1,:) = -model.t.^2.*exp(-model.t/0.02).*sin(2*pi*10*model.t); % Actividad a simular
% act(1,:) = act(1,:)/max(abs(act(1,:)));
% x = 10*(model.t-3);
% x = model.t*40-5;
% act(2,:) = exp(-(x).^2).*(-x.^5+x.^2);
% act(3,:) = -exp(-(model.t-0.15).^2/0.005).*sin(2*pi*5*model.t-3.5);
% act(1,:) = sin(2*pi*20*model.t+1.5);
% act(2,:) = sin(2*pi*40*model.t);
act(1,:) = 1;

[Laplacian, ~] = nip_neighbor_mat(model.cortex);
% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, 1*[15 20 15; 15 5 5; 15 -20 15; 15 5 25 ], ...
%         act,model.t);
[J, idx] = nip_simulate_activity(model.cortex,Laplacian, [25 0 15], ...
        act,model.t);
% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, size(act,1), ...
%         act,model.t);
fuzzy = nip_fuzzy_sources(model.cortex,1);
J = fuzzy*J;  
% rec_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
% subplot(1,2,1)
% plot(J(idx,:)')

model.y = model.L*J;
model.y = nip_addnoise(model.y,1000);

figure('Units','normalized','position',[0.1 0.1 0.3 0.3])
subplot(1,3,1)
nip_reconstruction3d(model.cortex,J,gca);

nbasis = model.Nd;
fuzzy = nip_fuzzy_sources(model.cortex,1);
basis = fuzzy(:,randi([1,model.Nd],nbasis,1));

A = model.L*basis;

[xx,status]=dalsql1(ones(nbasis,1), A, model.y, 0.000001);
% [xx,status]=dallrl1(zeros(128,1), 0, A, model.y, 0.5);

J_rec = basis*xx;
subplot(1,3,2)
nip_reconstruction3d(model.cortex,J_rec,gca)

Q = nip_lcmv(model.L,model.y);

[J_rec,~] = nip_loreta(model.y,model.L,diag(Q));
subplot(1,3,3)
nip_reconstruction3d(model.cortex,J_rec,gca)


