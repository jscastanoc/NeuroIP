% Script using sparse kalman

% Script for dictionaries generated using the principal components of the
% covariance in the time domain (London Style xD)
%% Init
clear; close all; clc;
nip_init();
verbose = true;

addpath('external');
addpath('external/kkmeans');
addpath('external/RWL1DF');
addpath('external/L1_homotopy');
setup_path;
% Load data, define sample rate, number of eeg channels, etc...
load_data;

% Simulate brain activity and generate pseudo EEG
gen_eeg;

% Show simulated activity
if verbose
    nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),[]);
    pause(0.01)
end

%% Preprocessing
% Depth bias compensation
% model.L = nip_depthcomp(model.L,0.1);

%% SVD freq decomposition

[y_proj, y_rec, Ur, Er] = nip_tempcomp(model.y, model.t, [0.1 40], 0.8);
Np = size(Ur,2);

if verbose
    figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
    Nh = ceil(sqrt(Np));
    Nw = ceil(sqrt(Np));
    ha = tight_subplot(Nh, Nw,  0.05, 0.1, 0.1);
    for i = 1:Np;
        axes(ha(i));
        plot(model.t,(y_proj(:,i)*Ur(:,i)')');
        title(strcat('Comp ',num2str(i)));
    end
    axes(ha(i+1))
    plot(model.t,model.y')
    title('Original')
    pause(0.01)
end

% Spatial dictionary in case we solve with, for example, S-FLEX
sp_dict = nip_fuzzy_sources(model.cortex,2,struct('dataset','montreal','save',true));
idx = find(sp_dict < 0.05*max(abs(sp_dict(:))));
sp_dict(idx) = 0;
sp_dict = sparse(sp_dict);

if verbose
    figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
    Nh = 1;
    Nw = ceil(sqrt(Np));
    ha = tight_subplot(Nh, Nw,  0.05, 0.1, 0.1);
end
%
aff = nip_fuzzy_sources(model.cortex,10,struct('dataset','montreal','save',true));
finalD =[];
for i = 1:Np
    D(:,i) = nip_sflex(y_proj(:,i),model.L,sp_dict, 1e-6);
    
    % Energy Threshold
    De = nip_energy(D(:,i));
    idx = find(De < 0.1*max(De));
    De(idx) = 0;
    for j = 1:3
        idx3 = (idx-1)*3+j;
        D(idx3,i) = 0;
    end
    D(:,i) = D(:,i)./norm(D(:,i));
    
    
    if verbose
        axes(ha(i));
        nip_reconstruction3d(model.cortex,D(:,i), []);
        pause(0.01)
    end
    
    
    % Clustering
    idx = find(De);
    affcr = aff(idx,idx);
    S = svd(affcr);
    
    Nclusters = numel(find(S >= 0.6*max(S(:))));
    label = knkmeans(affcr,Nclusters);
    
    
    labelfull = zeros(model.Nd/3,1);
    labelfull(idx) = label;
    
    finalDi = zeros(size(D(:,i),1),Nclusters+1);
    
    n = 1;
    for j = unique(labelfull)'
        idx1 = find(labelfull == j);
        for k = 1:3
            idx3 = (idx1-1)*3+k;
            finalDi(idx3,n) = D(idx3,i);
        end
        n = n+1;
    end
    Eaux(i) = Nclusters;
    finalD = [finalD,finalDi(:,2:end)];
end
% load('matlab.mat')
% finalD = [finalD,ones(model.Nd,1)];

for i=1:size(finalD,2)
    finalD(:,i) = finalD(:,i)./norm(finalD(:,i));
end


% idx = 1:3:model.Nd;
% n = 1;
% clear D
% D = zeros(model.Nd,size(finalD,2)*3);
% for i = 1:size(finalD,2);
%     fDe(:,i) = nip_energy(finalD(:,i));
%     for j = 0:2
%             D(idx+j,n) = fDe(:,i);
%             n = n+1;
%     end
% end

% Matrix with the perfectly located dictionaries mapped in to the lead
% field
Laux = model.L*finalD;


%% Compute hyperparameters and show results
fprintf('Computing Sparse Kalman...\n')
ha = nip_kalman_hyper_sparse(y_rec,Laux);

if verbose
    figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
    plot(model.t,ha')
    title('Temporal evolution of the hyperparameters')
    pause(0.01)
end

%% Final solution
clear J_est h;
% for i = 1:size(ha,1)
%     ha(i,:) = ha(i,:)/norm(ha(i,:));
% end

% for i = 1:size(Ur,2)
% %     Ecur = []
% %     for j = 1:length(Er)
% %         Ecur = [Ecur; repmat(Er(j),Eaux(j),1)]
% %     end
% %     Ecur = [Ecur; min(Ecur)];
% %     h =  (diag(1./Ecur))*ha;
%     hp = abs(ha*Ur(:,i))
%     Q = [];
%     for j = 1:size(D,2)
%         Q = [Q, hp(j)*finalD(:,j)];
%     end
%     Q = sum(Q,2);
%     
%     [J_est(:,i),~] = nip_loreta(y_proj(:,i),model.L,diag(Q));
% end


% J_est = J_est*Ur';

% for i = 1:model.Nt
%     i
%     [J_est(:,i),~] = nip_loreta(y_rec(:,i),model.L,diag(finalD*ha(:,i)));
% end