% Script PARAFAC
clear; close all; clc;
nip_init();
ltfatstart;
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
cfg.t = 0:1/cfg.fs:0.250; % Vector de tiempo
cfg.elec = elec;

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

act(1,:) = -model.t.^2.*exp(-model.t/0.02).*sin(2*pi*10*model.t); % Actividad a simular
act(1,:) = act(1,:)/max(abs(act(1,:)));
x = 10*(model.t-3);
x = model.t*40-5;
act(2,:) = exp(-(x).^2).*(-x.^5+x.^2);
% act(3,:) = -exp(-(model.t-0.15).^2/0.005).*sin(2*pi*5*model.t-3.5);
% act(4,:) = sin(2*pi*20*model.t+1.5);


[Laplacian, ~] = nip_neighbor_mat(model.cortex);
% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, 1*[15 20 15; 15 5 5; 15 -20 15; 15 5 25 ], ...
%         act,model.t);
[J, idx] = nip_simulate_activity(model.cortex,Laplacian, [5 0 15; 5 -5 15], ...
        act,model.t);
% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, size(act,1), ...
%         act,model.t);
fuzzy = nip_fuzzy_sources(model.cortex,1.5);
J = fuzzy*J;  
rec_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
subplot(1,2,1)
plot(J(idx,:)')

model.y = model.L*J;
model.y = nip_addnoise(model.y,10);

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
M = 500; % Frequency Res
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
Nfac =2; % Number of factors to decompose
[Factors,it,err,concord] = parafac(abs(c),Nfac,Opt,const);
figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
title('Factors')
subplot(1,3,1)
plot(Factors{1})
subplot(1,3,2)
plot(Factors{2})
subplot(1,3,3)
plot(Factors{3})


figure('Units','Normalized','Position',[0.3 0.1 0.3 0.5])

%%% Apply "Binary" mask on the original TF to extract each of the principal
%%% factors, then, reconstruct the time series of the signal directly from
%%% the TF.
for i = 1:Nfac
    for j=1:3
        Fact_rec{j} = Factors{j}(:,i);
    end
    S_rec = nmodel(Fact_rec,[],2);
    subplot(3,2,i)
    im_S(:,:) = S_rec(ch,:,:);
    imagesc(im_S)
    axis xy;
    subplot(3,2,i+2)
    idx_S = find(abs(S_rec)>0.1*max(max(max(abs(S_rec)))));
    bin_S = zeros(size(S_rec));
    bin_S(idx_S) = 1;
    c_inv = permute(c.*bin_S,[2 3 1]); % Apply "filtering" mask and reorder coordinates to performe the iSTFT
    f{i}=idgtreal(c_inv,'gauss',a,M)'; % Inverse STFT
    plot(f{i}(:,1:model.Nt)')
    
    subplot(3,2,i+4)
    timelock.individual =f{i}(:,1:model.Nt).^2;
    timelock.avg = f{i}(:,1:model.Nt).^2;
    ft_topoplotER(cfg, timelock); 
end
% Next: Find a solution with each of the f's

figure('Units','Normalized','Position',[0.3 0.1 0.5 0.3])
% Reconstruccion datos originales
Q = nip_lcmv(model.y,model.L); % Obtain beamformer with original data
Q = diag(Q);
[J_rec,~] = nip_loreta(model.y,model.L,Q);
subplot(1,3,1)
nip_reconstruction3d(model.cortex,sum(J_rec.^2,2),gca);

% Reconstruccion con cada modo identificado con PARAFAC
% for i = 1:Nfac
%     Q = nip_lcmv(f{i}(:,1:model.Nt),model.L);
%     Q = diag(Q);
%    [J_recSep{i},~] = nip_loreta(f{i}(:,1:model.Nt),model.L,Q);
%    subplot(1,3,i+1)
%    nip_reconstruction3d(model.cortex,sum(J_recSep{i}.^2,2),gca);
% end
