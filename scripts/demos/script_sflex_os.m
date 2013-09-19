clc, close all, clear;
% Script to solve the inverse problem using SFLEX. This particular script
% uses a modified version of a function written by Stefan Haufe (TU
% Berlin), the modification consist on the creation of the spatial
% dictionary outside the function (the original function creates the
% spatial dictionary inside of it)
%% Load Forward problem data
rng('default')
rng('shuffle');
nip_init();
warning off

load clab_full % Labels for 118 channels
clab_full = labels;
load clab_10_10; % Labels for 59 channels under 10 10 protocol
clab = clab_10_10;

load('data/sa_montreal.mat')
% data_name = 'icbm152b_sym';
% data_name = 'montreal';
% sa = prepare_sourceanalysis(clab, data_name);

temp = sa.V_cortex_coarse;
L = nip_translf(temp); % Leadfield matrix
L = L(find(ismember(clab_full,clab)),:);
clear temp 

cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = 200; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:0.5; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa cfg;

%% Simulation of neural activity
act = sin(2*pi*10*model.t); % Actividad a simular

% Simular actividad "act" en el dipolo más cercano a las coordenadas [30
% -20 30]
[J, active] = nip_simulate_activity(model.cortex.vc,[30 -20 30], act, randn(1,3), model.t);

% Hay dispersion de la actividad, entre mas grande el numero, mas dispersa es la actividad    
fuzzy = nip_fuzzy_sources(model.cortex,0.001);
index = (1:3:model.Nd);

for i = 0:2
    J(index+i,:) = fuzzy*J(index+i,:); % J simulado FINAL
end

% Obtener el eeg correspondiente a la simulacion
clean_y = model.L*J;

% Anadir ruido
snr = 10;
model.y = nip_addnoise(clean_y, snr);

%% Pre-process 
model.L = nip_depthcomp(model.L,0.2); % Depth bias compensation for the lead field matrix

% Construction of the spatial basis functions
nbasis = size(model.L,2)/3; % Number of basis functions
iter_basis = [1.5]; % Width of the spatial basis functions
basis = [];
for i = iter_basis
    fuzzy = nip_fuzzy_sources(model.cortex,i);
    basisn = fuzzy(:,randi([1,model.Nd/3],nbasis,1));
    basisn = basisn/norm(basisn(:),1);
    basis = [basis basisn];
end

%% Solution of the inverse problem
% model.L = nip_translf(model.L);
index = (1:3:model.Nd);
for i = 0:2
    L(:,index+i) = model.L(:,index+i)*basis; % J simulado FINAL
end
% L = nip_translf(L);
[xx,~] = sflex_cortical_dal(model.y,L,[],struct('eps',0.25,'B',eye(model.Nd)));
% for i = 0:2
%     xx(index+i,:) = basis*xx(index+i,:); % J simulado FINAL
% end

% J_rec = nip_translf(permute(xx,[3 1,2]))';
J_rec = permute(xx,[1 3,2]);
index = (1:3:model.Nd);
for i = 0:2
    J_rec(index+i,:) = basis*J_rec(index+i,:); % J simulado FINAL
end
% J_rec = nip_trans_solution(J_rec);
% [J_rec,~] = nip_sflex(model.y,model.L,basis);
%%%%%%%%%%%%%%%%%
% Visualization %
%%%%%%%%%%%%%%%%%
% Simulation - Temporal
figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
plot(model.t,J')
xlabel('Time')
ylabel('Amplitude')

% Simulation - Spatial
figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
nip_reconstruction3d(model.cortex, sqrt(sum(J.^2,2)), struct('axes',gca));
hold on
scatter3(model.cortex.vc(active,1),model.cortex.vc(active,2),model.cortex.vc(active,3),'ok','filled');


%%% Reconstruction - Temporal
figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
plot(model.t,J_rec')
xlabel('Time')
ylabel('Amplitude')

% Reconstruction - Spatial
figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
nip_reconstruction3d(model.cortex, sqrt(sum(J_rec.^2,2)), struct('axes',gca));
hold on
h = scatter3(model.cortex.vc(active,1),model.cortex.vc(active,2),model.cortex.vc(active,3),'ok','filled');
