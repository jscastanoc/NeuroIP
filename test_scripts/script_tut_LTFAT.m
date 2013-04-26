% Script LTFAT
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
cfg.fs = 500; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:0.250; % Vector de tiempo

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

act(1,:) = -model.t.^2.*exp(-model.t/0.02).*sin(2*pi*10*model.t); % Actividad a simular
act(1,:) = act(1,:)/max(abs(act(1,:)));
x = 20*(model.t-3);
x = model.t*80-5;
act(2,:) = exp(-(x).^2).*(-x.^5+x.^2);
act(3,:) = -exp(-(model.t-0.15).^2/0.005).*sin(2*pi*5*model.t-3.5);
% act(4,:) = sin(2*pi*10*model.t);

actorg = sum(act,1);
actfull = nip_addnoise(actorg,20);
a = 1;
M = 126 ;

c = dgtreal(actfull,'gauss',a,M);
figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
plotdgtreal(c,a,M,'linsq');
disp('Number of time shifts in transform:')
N = size(c,2)

disp('Length of transform:')
L = N*a


figure('Units','Normalized','Position',[0.1 0.1 0.3 0.3])
subplot(2,1,1)
plot(actfull)
hold on
plot(actorg,'r')

hold
subplot(2,1,2)
c_sgram=sgram(actfull,model.fs,'lin');