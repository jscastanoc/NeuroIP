% Script for dictionaries generated using the principal components of the
% covariance in the time domain (London Style xD)
%% Init
clear; close all; clc;

verbose = true;

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
model.L = nip_depthcomp(model.L,0.1);

%% SSA decomposition

nf = 30;
for i = 1:model.Nc
    y_decomp(i,:,:) = MSSA_Hankel(model.y(i,:),nf,nf);
end

break
Np = 2;

ybase = zeros(1,model.Nt);
for i = 1:Np
    ybase(:,:) = mean(y_decomp(:,:,i));
    y_proj(:,i+1) = model.y*ybase';
end

break
% if verbose
%     figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
%     Nh = ceil(sqrt(Np));
%     Nw = ceil(sqrt(Np));
%     ha = tight_subplot(Nh, Nw,  0.05, 0.1, 0.1);
%     for i = 1:Np;
%         axes(ha(i));
%         plot(model.t,y_proj(:,:,i));
%         title(strcat('Comp ',num2str(i)));
%     end
%     axes(ha(i+1))
%     plot(model.t,model.y')
%     title('Original')
%     pause(0.01)
% end

% Spatial dictionary in case we solve with, for example, S-FLEX
sp_dict = nip_fuzzy_sources(model.cortex,2);
idx = find(sp_dict < 0.05*max(abs(sp_dict(:))));
sp_dict(idx) = 0;
sp_dict = sparse(sp_dict);
for i = 1:Np
    %     D(:,i) = nip_lcmv(y_proj(:,i),model.L);
    D(:,i) = nip_sflex(y_proj(:,i),model.L,sp_dict, 1e-6);
    %     [aux,~] = nip_tfmxne_port(Ye(:,:,i),model.L,[])
    %     D(:,i) = sqrt(sum(aux.^2,2));
%     D(:,i) = D(:,i) - mean(D(:,i));
    D(:,i) = D(:,i)./norm(D(:,i));
    
    if verbose
        axes(ha(i));
        nip_reconstruction3d(model.cortex,D(:,i), []);
        pause(0.01)
    end
end

L = model.L*D;

%% Compute hyperparameters and show results
h = nip_kalman_hyper(model.y,L);
if verbose
    figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
    plot(model.t,h')
    title('Temporal evolution of the hyperparameters')
    pause(0.01)
end

%% solve