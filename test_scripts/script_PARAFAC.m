% Script PARAFAC
clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 4000; % 1000, 2000, 4000, 8000

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('data/montreal',num2str(Nd),'_full3shell.mat'))

cfg.L = L;
cfg.cortex = cortex_mesh;
cfg.fs = 1000; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:0.250; % Vector de tiempo

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

act(1,:) = -model.t.^2.*exp(-model.t/0.02).*sin(2*pi*10*model.t); % Actividad a simular
act(1,:) = act(1,:)/max(abs(act(1,:)));
x = 10*(model.t-3);
x = model.t*40-5;
act(2,:) = exp(-(x).^2).*(-x.^5+x.^2);
act(3,:) = -exp(-(model.t-0.15).^2/0.005).*sin(2*pi*5*model.t-3.5);
act(4,:) = sin(2*pi*6*model.t);


[Laplacian, ~] = nip_neighbor_mat(model.cortex);
[J, idx] = nip_simulate_activity(model.cortex,Laplacian, 1*[15 20 15; 15 5 5; 15 -20 15; 15 5 25 ], ...
        act,model.t);
fuzzy = nip_fuzzy_sources(model.cortex,1);
J = fuzzy*J;  
rec_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
subplot(1,2,1)
plot(J(idx,:)')

model.y = model.L*J;
model.y = nip_addnoise(model.y,10);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time Frequency decomposition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fb = 1.5;
f = 0.5:1:30; % Bins of the time frequency map
Fs = model.fs;
Nf = length(f);
for j=1:Nf
    f0 =f(j);
    sigma_f = f0/7;
    sigma_t = 1/(2*pi*sigma_f);
    t = -2*sigma_t:1/Fs:2*sigma_t;
    A = (sigma_t)^(-1/2)*pi^(-1/4);
    w = A.*exp(2*1i*pi*f0*t).*exp(-(t.^2)/(2*sigma_t^2));
    parfor i = 1:model.Nc
        S(i,j,:) = abs(conv(model.y(i,:),w,'same')); %% Time frequency maps
    end
end

%%%%%%%%%%%
% Show TF %
%%%%%%%%%%%
figure('Units','normalized','position',[0.1 0.1 0.3 0.5]);

for i = 1:model.Nc
    img_ch(:,:) = S(i,:,:);
    imagesc(img_ch,[min(min(min(S))),max(max(max(S)))])
    axis xy;
    pause(1)
end

% Parafac Decomposition
Opt(1) = 1e-6; Opt(2) = 1; Opt(3) = 0; Opt(4) = 0; Opt(5) = 10; Opt(6) = 2500;
const = [2,2,2];
Nfac =4; % Number of factors to decompose
[Factors,it,err,concord] = parafac(S,Nfac,Opt,const);


S_rec = nmodel(Factors,[],2);



