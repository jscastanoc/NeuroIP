% Testbench Dynamic-Projected
clear; close all; clc;
nip_init();
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 2000; % 1000, 2000, 4000, 8000

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('../../data/nocons_or/montreal',num2str(Nd),'_full1shell.mat'))

cfg.L = lf;
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

for i = 1:3
    acti{i} = act; % Actividad a simular
end
act = acti;
clear acti;

[~,index_t1] = max(act{1}(1,:));
[~,index_t2] = max(act{1}(2,:));

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

% SNR = [-10 0 10];
% NDIP = [1,2,3];
SNR = [0];
NDIP = [2];
NITER = 1;
dir = '/home/jscastanoc/Desktop/Results_TV_Prios/T3';
save_res = false;

graphic_mode = true;

if graphic_mode
    NITER = 1;
    save_res = false;
    fig_sp = figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
end


for iter = 1:NITER
    for ndip = NDIP
        [J, idx] = nip_simulate_activity(model.cortex,  ndip, ...
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
            methods = {'S-FLEX'}
            [yp, yr, Ur, ~] = nip_tempcomp(model.y, model.t, [0 60], 0.7);
            Nr = size(Ur,2);
            for m_idx = 1:length(methods);
                method = methods{m_idx};
                tic;
                switch method
                    case 'BMF'
                        patch_set = BMFUr(yp,model.L);
                    case 'ARD'                        
                        patch_set = MSPUr(yp, model.L, basis, method);
                    case 'GS'
                        patch_set = MSPUr(yp, model.L, basis, method);
                    case 'S-FLEX'
                        patch_set = SFLEXUr(yp, model.L,basis);            
                end
                for k = 1:Nr
                    h(:,k) = nip_spm_priors(yp(:,k), model.L, patch_set, Qe, 'ARD');
                    Q = zeros(model.Nd, 1);
                    for i = 1:size(patch_set,2)
                        Q = Q + h(i,k)*patch_set(:,i);
                    end
                    Q = diag(Q);
                    [J_estd(:,k),~] = nip_spm_inversion(yp(:,k),model.L,Q,Qe);
                end
                J_est = J_estd*Ur';
                et = toc;
                
                if save_res
                    file_name = strcat(dir,method,'_',num2str(snr),'_',num2str(ndip),...
                        '_',num2str(iter));
                    y = model.y;
                    save(file_name,'J','J_est','y','et','h','patch_set','Ur');
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