% Testbench Static

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
act(3,:) = -exp(-(model.t-0.15).^2/0.005).*sin(2*pi*5*model.t-3.5);
act(4,:) = sin(2*pi*20*model.t+1.5);


[Laplacian, ~] = nip_neighbor_mat(model.cortex);
% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, 1*[15 20 15; 15 5 5; 15 -20 15; 15 5 25 ], ...
%         act,model.t);

% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, size(act,1), ...
%         act,model.t);


% Spatial Dictionary for ARD, GS and S-FLEX
nbasis = 500;
fuzzy = nip_fuzzy_sources(model.cortex,0.02);
basis1 = fuzzy(:,randi([1,model.Nd],nbasis,1));
fuzzy = nip_fuzzy_sources(model.cortex,0.027);
basis2 = fuzzy(:,randi([1,model.Nd],nbasis,1));
fuzzy = nip_fuzzy_sources(model.cortex,0.035);
basis3 = fuzzy(:,randi([1,model.Nd],nbasis,1));

basis = [basis1, basis2, basis3];


fuzzy = nip_fuzzy_sources(model.cortex,0.02);


Qe = eye(model.Nc);

SNR = [-10 0 10];
NDIP = [1,2,3];
% SNR = [-10];
% NDIP = [1];
NITER = 30;
dir = '/home/jscastanoc/Desktop/Results_TV_Prios/T2';
save_res = true;

graphic_mode = false;
if graphic_mode
    NITER = 1;
    save_res = false;
    fig_sp = figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
end


for iter = 1:NITER
    for ndip = NDIP
        [J, idx] = nip_simulate_activity(model.cortex, Laplacian,  ndip, ...
            act(1:ndip,:),model.t);
        J = fuzzy*J;
        J = sparse(J);
        y_clean = model.L*J;
        
        if graphic_mode
           figure(fig_sp)
           subplot(2,3,1)
           nip_reconstruction3d(model.cortex, sqrt(sum(J.^2,2)),gca)
        end
        for snr = SNR
            model.y = nip_addnoise(y_clean,snr);
%             methods = {'BMF','ARD','GS','S-FLEX'};
            methods = {'S-FLEX','GS'}
            [yp, yr, Ur, ~] = nip_tempcomp(model.y, model.t, [0 60], 0.7);
            Nr = size(Ur,2);
            for m_idx = 1:length(methods);
                method = methods{m_idx};
                tic;
                switch method
                    case 'BMF'
                        BMFpatch = BMFUr(yp,model.L);
                        h = nip_spm_priors(yp, model.L, BMFpatch, Qe, 'GS');
                        h = h/max(h);
                        Q = zeros(model.Nd, 1);
                        for i = 1:size(BMFpatch,2)
                            Q = Q + h(i)*BMFpatch(:,i);
                        end
                        Q = diag(Q);
                    case 'ARD'                        
                        MSPpatch = MSPUr(yp, model.L, basis,method);
                        h = nip_spm_priors(yp, model.L, MSPpatch, Qe, 'GS');
                        h = h/max(h);
                        Q = zeros(model.Nd, 1);
                        for i = 1:size(MSPpatch,2)
                            Q = Q + h(i)*MSPpatch(:,i);
                        end
                        Q = diag(Q);
                    case 'GS'
                        MSPpatch = MSPUr(yp, model.L, basis,method);
                        h = nip_spm_priors(yp, model.L, MSPpatch, Qe, 'GS');
                        h = h/max(h);
                        Q = zeros(model.Nd, 1);
                        for i = 1:size(MSPpatch,2)
                            Q = Q + h(i)*MSPpatch(:,i);
                        end
                        Q = diag(Q);
                    case 'S-FLEX'
                        SFLEXpatch = SFLEXUr(yp, model.L,basis);
                        h = nip_spm_priors(yp, model.L, SFLEXpatch, Qe, 'GS');
                        h = h/max(h);
                        Q = zeros(model.Nd, 1);
                        for i = 1:size(SFLEXpatch,2)
                            Q = Q + h(i)*SFLEXpatch(:,i);
                        end
                        Q = diag(Q);
%                     	[J_est,~] = nip_loreta(yp,model.L,speye(model.Nd));
            
                end
                if sum(strcmp(method,{'BMF','ARD','GS','S-FLEX'}))
                    [J_est,~] = nip_spm_inversion(yp,model.L,Q,Qe);
                end
                J_est = J_est*Ur';
                et = toc;
                
                if save_res
                    file_name = strcat(dir,method,'_',num2str(snr),'_',num2str(ndip),...
                    '_',num2str(iter));
                    y = model.y;
                    save(file_name,'J','J_est','y','et');
                end
                if graphic_mode
                    subplot(2,3,m_idx+1)
                    nip_reconstruction3d(model.cortex, sqrt(sum(J_est.^2,2)),gca)
                    title(method)
                end
                clear Q sf_patch bmf_patch;
            end
        end
    end        
end