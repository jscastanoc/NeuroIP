% Script BASIC
clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('data/sa_montreal.mat'))

cfg.L = nip_translf(sa.V_cortex_coarse);
cfg.cortex = sa.cortex_coarse;
cfg.fs = 200; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:1; % Vector de tiempo

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec

% Simulated activity
act = sin(2*pi*10*model.t);

% Simular actividad "act" en el dipolo más cercano a las coordenadas [30
% -20 30]
[J, ~] = nip_simulate_activity(model.cortex.vc,[30 -20 30], act, randn(1,3), model.t);

% Hay dispersion de la actividad, entre mas grande el numero, mas dispersa es la actividad
fuzzy = nip_fuzzy_sources(model.cortex,1);
index = (1:3:model.Nd);

for i = 0:2
    J(index+i,:) = fuzzy*J(index+i,:); % J simulado FINAL
end

% Obtener el eeg correspondiente a la simulacion
clean_y = model.L*J;
snr = 10;
model.y = nip_addnoise(clean_y, snr);



% Depth compensation
depth = 'Lnorm';
switch depth
    case 'none'
        L = model.L;
    case 'Lnorm'        
        gamma = 0.9;
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;        
end
clear extras;

%% SOLUTION %%

method = 'STOUT';
switch method
    case 'LORETA'
        Q = eye(model.Nd); %Matriz de covarianza apriori
        [J_est, extras] = nip_loreta(model.y, L, Q);
    case 'TFMxNE'
        options.spatial_reg = 80;
        options.temp_reg =  0.5;
        options.iter = 50;
        options.tol = 2e-2;
        [J_est, extras] = nip_tfmxne_port(model.y, L, options);
    case 'STOUT' 
        % Spatial dictionary
        sigma = 1;
        B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
        
        % Options for the inversion
        options.spatial_reg = 80;
        options.temp_reg =  0.5;
        options.iter = 50;
        options.tol = 2e-2;
        [J_est, extras] = nip_stout(model.y, L, B, options);
    case 'SFLEX'
        % Spatial dictionary
        sigma = 1;        
        B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
        
        % Options for the inversion
        reg_par = 1;
        [J_est, extras] = nip_sflex(model.y, L, B, reg_par);
end


if ~strcmp(depth,'none')
    J_est = nip_translf(J_est');
    for i = 1:3
        J_est(:,:,i) = (Winv(:,:,i)*J_est(:,:,i)')';
    end
    J_est = nip_translf(J_est)';
end

%% ERROR MEASURES %% 
r = 5;
[out, Ms, Mr] = nip_error_sai(model.cortex, J,J_est, r);
out
       
dist = nip_fuzzy_sources(model.cortex,[], ...
    struct('save',1,'dataset','montreal','calc','dist'));


emd = nip_emd(J,J_est,dist);

%% Visualization %%
%%% Visualizacion de la simulacion %%%
%Temporal
figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
plot(model.t,J')
xlabel('Time')
ylabel('Amplitude')


figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
nip_reconstruction3d(model.cortex, sqrt(sum(J.^2,2)), struct('axes',gca)); hold on;
idx_localmax = find(Ms);
scatter3(model.cortex.vc(idx_localmax,1),...
    model.cortex.vc(idx_localmax,2),...
    model.cortex.vc(idx_localmax,3),...
    'filled','MarkerFaceColor','k')
%%% Visualizacion de la reconstruccion %%%
%Temporal
figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
plot(model.t,J_est')
xlabel('Time')
ylabel('Amplitude')


figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
nip_reconstruction3d(model.cortex, sqrt(sum(J_est.^2,2)),  struct('axes',gca)); hold on;
idx_localmax = find(Mr);
scatter3(model.cortex.vc(idx_localmax,1),model.cortex.vc(idx_localmax,2),model.cortex.vc(idx_localmax,3),'filled')

