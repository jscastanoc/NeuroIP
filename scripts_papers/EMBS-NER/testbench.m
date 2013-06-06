% Testbench NER

clear; close all; clc;
nip_init();
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

% Spatial Dictionary for ARD and S-FLEX
nbasis = 1024;
fuzzy = nip_fuzzy_sources(model.cortex,1);
fuzzy_patch = nip_fuzzy_sources(model.cortex,1.5);
basis = fuzzy_patch(:,randi([1,model.Nd],nbasis,1));

Qe = eye(model.Nc);

SNR = [-10 0 10];
NDIP = [1,2,3,4];
% SNR = [-10];
% NDIP = [1];
NITER = 30;
dir = '/home/jscastanoc/Desktop/Results_NER_1010/';
for iter = 1:NITER
    for ndip = NDIP
        [J, idx] = nip_simulate_activity(model.cortex, Laplacian, ndip, ...
            act(1:ndip,:),model.t);
        J = fuzzy*J;
        J = sparse(J);
        y_clean = model.L*J;
        
        for snr = SNR
            model.y = nip_addnoise(y_clean,snr);
            methods = {'BMF','BMF-TF','ARD','ARD-TF'};
%             methods = {'ARD'}
            for m_idx = 1:length(methods);
                method = methods{m_idx};
                tic;
                switch method
                    case 'BMF'
                        Q = diag(nip_lcmv(model.y, model.L));
                    case 'BMF-TF'
                        bmf_patch = BMFTF(model.y,model.L);
                        h = nip_spm_solvers(model.y, model.L, bmf_patch, Qe, 'ARD');
                        h = h/max(h);
                        Q = zeros(model.Nd, 1);
                        for i = 1:size(bmf_patch,2)
                            Q = Q + h(i)*bmf_patch(:,i);
                        end
                        Q = diag(Q);
                    case 'ARD'
                        h = nip_spm_solvers(model.y, model.L, basis, Qe, method);
                        h = h/max(h);
                        Q = zeros(model.Nd, 1);
                        for i = 1:size(basis,2)
                            Q = Q + h(i)*basis(:,i);
                        end
                        Q = diag(Q);
                    case 'S-FLEX-TF'
                        sf_patch = SFLEXTF(model.y,model.L,basis);
                        h = nip_spm_solvers(model.y, model.L, sf_patch, Qe, 'ARD');
                        h = h/max(h);
                        Q = zeros(model.Nd, 1);
                        for i = 1:size(sf_patch,2)
                            Q = Q + h(i)*sf_patch(:,i);
                        end
                        Q = diag(Q);
                    case 'ARD-TF'
                        ard_patch = ARDTF(model.y,model.L,basis)
                        h = nip_spm_solvers(model.y, model.L, ard_patch, Qe, 'ARD');
                        h = h/max(h);
                        Q = zeros(model.Nd, 1);
                        for i = 1:size(ard_patch,2)
                            Q = Q + h(i)*ard_patch(:,i);
                        end
                        Q = diag(Q);
                end
                [J_est,~] = nip_loreta(model.y,model.L,Q);
                et = toc;
                file_name = strcat(dir,method,'_',num2str(snr),'_',num2str(ndip),...
                    '_',num2str(iter));
                y = model.y;
                save(file_name,'J','J_est','y','et');
                clear Q sf_patch bmf_patch;
            end
        end
    end        
end