% Script SFLEX

clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 2000; % 1000, 2000, 4000, 8000

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('../data/nocons_or/montreal',num2str(Nd),'_full1shell.mat'))

cfg.L = lf;
cfg.cortex = cortex_mesh;
cfg.fs = 200; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:0.25; % Vector de tiempo
cfg.elec = elec;

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

act = 100*sin(2*pi*10*model.t); % Actividad a simular



% [J, idx] = nip_simulate_activity(model.cortex, 1*[15 20 15; 15 5 5; 15 -20 15; 15 5 25 ], ...
%         act,model.t);
[J, idx] = nip_simulate_activity(model.cortex, [25 0 15], ...
        act,randn(size(act,1),3),model.t);
% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, size(act,1), ...
%         act,model.t);
fuzzy = nip_fuzzy_sources(model.cortex,0.11); 

index = (1:3:model.Nd);

for i = 0:2
    J(index+i,:) = fuzzy*J(index+i,:); % J simulado FINAL
end

% rec_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
% subplot(1,2,1)
% plot(J(idx,:)')

model.y = model.L*J;
model.y = nip_addnoise(model.y,10);

figure('Units','normalized','position',[0.1 0.1 0.3 0.3])
subplot(1,3,1)
nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),gca);


nbasis = 256;

iter_basis = [0.05 0.1 0.2];
basis = [];
for i = iter_basis
    fuzzy = nip_fuzzy_sources(model.cortex,i);
    basisn = fuzzy(:,randi([1,model.Nd/3],nbasis,1));
    basisn = basisn/norm(basisn(:),1);
    basis = [basis basisn ];
end



basis = reshape(repmat(basis(:)',3,1),model.Nd,nbasis*length(iter_basis));


[J_rec,~] = nip_sflex(model.y,model.L,basis,0.1);

subplot(1,3,2)
nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca)

Q = nip_lcmv(model.y,model.L);

[J_reclor,~] = nip_loreta(model.y,model.L,diag(Q));
subplot(1,3,3)
nip_reconstruction3d(model.cortex,sqrt(sum(J_reclor.^2,2)),gca)


