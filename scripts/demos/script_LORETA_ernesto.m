% Script BASIC
clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cargar datos(lead field, mesh del cerebro etc...
lpath = '/home/jscastanoc/Dropbox/Solo/Master/NeurocentroSoft/testbench_matlab/data/';
head = 'headN';
[cfg.L cfg.cortex] = nip_ernestoLF(lpath,head);
cfg.fs = 200; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:1; % Vector de tiempo

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

act = sin(2*pi*10*model.t); % Actividad a simular

% Simular actividad "act" en el dipolo más cercano a las coordenadas [30
% -20 30]
[J, ~] = nip_simulate_activity(model.cortex.vc, 1, act, randn(1,3), model.t);

% Hay dispersion de la actividad, entre mas grande el numero, mas dispersa es la actividad    
% fuzzy = nip_fuzzy_sources(model.cortex,0.1);
% index = (1:3:model.Nd);
% 
% for i = 0:2
%     J(index+i,:) = fuzzy*J(index+i,:); % J simulado FINAL
% end

% Obtener el eeg correspondiente a la simulacion
clean_y = model.L*J;

% Depth compensation
depth = 'Lnorm'; % it can be none, Lnorm or sLORETA-based depth compensation
switch depth
    case 'none'
        L = model.L;
        Winv = [];
    case 'Lnorm'
        gamma = 0.7;
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;
end

% Anadir ruido
snr = 10;
model.y = nip_addnoise(clean_y, snr);

% Average reference
transM = eye(model.Nc)-(1/model.Nc)*ones(model.Nc);
model.y = transM*model.y;
model.L = transM*model.L;


%%%%%%%%%%%%%%
% Estimacion %
%%%%%%%%%%%%%%
% Estimar actividad (en este caso LORETA por que se usa el laplaciano
% espacial para hallar la matriz de covarianza
Q = speye(model.Nd); %Matriz de covarianza apriori

%%% Beamformer as prior covariance
% Q = diag(nip_lcmv(model.y,model.L)); %requires spm_inv (SPM toolbox)

[J_est, extras] = nip_loreta(model.y, model.L, 'cov', Q, 'Winv', Winv);


%%%%%%%%%%%%%%%%%
% OTHER METHODS %
%%%%%%%%%%%%%%%%%
% sigma = 1.5; % Width of the spatial basis functions
% B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
% B = nip_blobnorm(B,'norm',2);
% [J_est, ~] = nip_sflex(y, L, B); % Needs DAL toolbox


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
nip_reconstruction3d(model.cortex, J(:,t_0), struct('axes',gca));

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
nip_reconstruction3d(model.cortex, J_est(:,t_0),  struct('axes',gca));