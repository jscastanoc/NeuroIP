% Script PARAFAC
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
cfg.t = 0:1/cfg.fs:0.250; % Vector de tiempo
cfg.elec = elec;

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

%%%%%%%%%%%%%%%%%%%%%%
% Actividad a simular
%%%%%%%%%%%%%%%%%%%%%%%%
act(1,:) = -model.t.^2.*exp(-model.t/0.02).*sin(2*pi*10*model.t); 
act(1,:) = act(1,:)/max(abs(act(1,:)));
x = 10*(model.t-3);
x = model.t*40-5;
act(2,:) = exp(-(x).^2).*(-x.^5+x.^2);
% act(3,:) = -exp(-(model.t-0.15).^2/0.005).*sin(2*pi*5*model.t-3.5);
% act(4,:) = sin(2*pi*20*model.t+1.5);


[Laplacian, ~] = nip_neighbor_mat(model.cortex);
% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, 1*[15 20 15; 15 5 5; 15 -20 15; 15 5 25 ], ...
%         act,model.t);
[J, idx] = nip_simulate_activity(model.cortex,Laplacian, [5 0 15; 5 -3 15], ...
        act,model.t);
% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, size(act,1), ...
%         act,model.t);
fuzzy = nip_fuzzy_sources(model.cortex,1);
J = fuzzy*J;  
rec_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
subplot(1,2,1)
plot(J(idx,:)')

model.y = model.L*J;
model.y = nip_addnoise(model.y,-10);

clear cfg;
cfg.layout = 'EEG1010.lay';
% cfg.projection = 'inverse';
% cfg.overlap = 'keep';
% cfg.rotate = 3;
% data.elec = model.elec.elecpos;
% data.label = model.elec.label;
lay = ft_prepare_layout(cfg);
clear cfg

ch = find(ismember(model.elec.label,lay.label));
timelock.dimord ='chan_time';
timelock.individual(1,:,:) = model.y.^2;
timelock.avg = model.y.^2;
timelock.label = model.elec.label(ch);
timelock.time = model.t;
cfg.layout = lay;
figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
subplot(1,2,1)
nip_reconstruction3d(model.cortex,sum(J.^2,2),gca)
% nip_reconstruction3d(model.cortex,J(:,50),gca)
subplot(1,2,2)
ft_topoplotER(cfg, timelock);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Frequency decomposition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 5; % Time shift
M = 200; % Frequency Res
c = dgtreal(model.y','gauss',a,M);
c = permute(c,[3 1 2]);

%%%%%%%%%%%
% Show TF %
%%%%%%%%%%%
ch = 30; %% Show the TF of this channel
figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
subplot(2,1,1)
plot(model.t,model.y(ch,:))
subplot(2,1,2)
% c_sgram=sgram(model.y(ch,:),model.fs,'lin','wlen',30);
cc(:,:) = c(ch,:,:);
plotdgtreal(cc,a,M,'linsq');



% Parafac Decomposition
Opt(1) = 1e-6; Opt(2) = 1; Opt(3) = 0; Opt(4) = 0; Opt(5) = 10; Opt(6) = 2500;
const = [2,2,2];
Nfac =3; % Number of factors to decompose
[Factors,it,err,concord] = parafac(abs(c),Nfac,Opt,const);
figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
title('Factors')
subplot(1,3,1)
plot(Factors{1})
subplot(1,3,2)
plot(Factors{2})
subplot(1,3,3)
plot(Factors{3})



% SFLEX Preprocessing
nbasis = model.Nd;
fuzzy = nip_fuzzy_sources(model.cortex,0.8);
basis = fuzzy(:,randi([1,model.Nd],nbasis,1));
A = model.L.^2*basis;


% Reconstruction with each factor
for i = 1:Nfac
    [xx,status]=dalsql1(ones(nbasis,1), A, Factors{1}(:,i), 0.001);
    J_rec{i}  = basis*xx;
    figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
    nip_reconstruction3d(model.cortex,J_rec{i},gca);
end

figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
nip_reconstruction3d(model.cortex,nip_lcmv(model.y,model.L)',gca);
