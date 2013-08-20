
clear; close all; clc;
nip_init();

addpath('../../external/source_toolbox/haufe/')
addpath('../../external/source_toolbox/nolte/')
addpath('../../external/source_toolbox/simulations/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 1000; % 1000, 2000, 4000, 8000

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('../../data/nocons_or/montreal',num2str(Nd),'_full1shell.mat'))

cfg.L = lf;
cfg.cortex = cortex_mesh;
cfg.fs = 200; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:0.25; % Vector de tiempo
cfg.elec = elec;

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

act = sin(2*pi*10*model.t); % Actividad a simular



% [J, idx] = nip_simulate_activity(model.cortex.vertices, 1*[15 20 15; 15 5 5; 15 -20 15; 15 5 25 ], ...
%         act,model.t);
[J, idx] = nip_simulate_activity(model.cortex.vertices, [25 0 15], ...
        act,randn(size(act,1),3),model.t);

fuzzy = nip_fuzzy_sources(model.cortex,3); 

index = (1:3:model.Nd);

for i = 0:2
    J(index+i,:) = fuzzy*J(index+i,:); % J simulado FINAL
end

rec_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);

model.y = model.L*J;
model.y = nip_addnoise(model.y,0);

rec_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
subplot(1,3,1)
nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),gca);

nbasis = 1000;
iter_basis = [2];
basis = [];
n = 1;
group = [];
for i = iter_basis
    fuzzy = nip_fuzzy_sources(model.cortex,i);
    basisn = fuzzy(:,randi([1,model.Nd/3],nbasis,1));
    basisn = basisn/norm(basisn(:),1);
    basis = [basis basisn ];
    group = [group n*ones(1,nbasis)];
    n = n+1;
end
basis = nip_blobnorm(basis,group,struct('norm',1,'norm_group',true));

model.L = nip_depthcomp(model.L,0.2);

% [P, out] = loreta(model.y, nip_translf(model.L), model.cortex.vertices,[]);
% J_rec = permute(P(:,:,:,40),[3 1 2]);
% J_rec = nip_translf(J_rec)';

options.iter = 50;
options.spatial_reg = 2;
options.temp_reg = 1;
options.tol = 1e-2;
[J_rec,~] = nip_sflex_tfmxne(model.y,model.L,basis,options);

figure;
plot(J_rec')
title('TF')

figure(rec_fig)
subplot(1,3,2)
nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);

dist = nip_emd(sqrt(sum(J.^2,2)),sqrt(sum(J_rec.^2,2)),model.cortex);
fprintf('EMD ST spatial: %.4d \n', dist)
dist = nip_emd(sqrt(sum(J.^2,2)),sqrt(sum(J_rec.*J,2)),model.cortex);
fprintf('EMD ST temp: %.4d \n', dist)

Q = nip_lcmv(model.y,model.L);
Q = Q/max(Q);

[J_reclor,~] = nip_loreta(model.y,model.L,diag(Q));
subplot(1,3,3)
nip_reconstruction3d(model.cortex,sqrt(sum(J_reclor.^2,2)),gca);
figure;
plot(J_reclor')
title('LORETA')

dist = nip_emd(sqrt(sum(J.^2,2)),sqrt(sum(J_reclor.^2,2)),model.cortex);
fprintf('EMD LOR spatial: %.4d \n', dist)
dist = nip_emd(sqrt(sum(J.^2,2)),sqrt(sum(J_reclor.*J,2)),model.cortex);
fprintf('EMD LOR temp: %.4d \n', dist)
