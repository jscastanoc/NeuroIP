% Script for dictionaries generated using the principal components of the
% covariance in the time domain (London Style xD)
%% Init
% clear; close all; clc;
close all;
nip_init();
verbose = true;

% Load data, define sample rate, number of eeg channels, etc...
load_data;

% Simulate brain activity and generate pseudo EEG
gen_eeg;

opt3d = struct('view',[90 0]);
% Show simulated activity
if verbose
    nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),opt3d);
    pause(0.01)
end
break
%% Preprocessing
% Depth bias compensation
% model.L = nip_depthcomp(model.L,0.1);

%% SVD freq decomposition

[y_proj, y_rec, Ur, Er] = nip_tempcomp(model.y, model.t, [0.1 40], 2);
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
savefig('eeg',gcf,'eps')

% Spatial dictionary in case we solve with, for example, S-FLEX
sp_dict = nip_fuzzy_sources(model.cortex,2,struct('dataset','montreal','save',true));
idx = find(sp_dict < 0.05*max(abs(sp_dict(:))));
sp_dict(idx) = 0;
sp_dict = sparse(sp_dict);



aff = nip_fuzzy_sources(model.cortex,4.5,struct('dataset','montreal','save',true));
finalD =[];
for i = 1:Np
    D(:,i) = nip_lcmv(y_proj(:,i),model.L);
        D(:,i) = nip_sflex(y_proj(:,i),model.L,sp_dict, 1e-6);
    %     [aux,~] = nip_tfmxne_port(Ye(:,:,i),model.L,[])
    %     D(:,i) = sqrt(sum(aux.^2,2));
    %     D(:,i) = D(:,i) - mean(D(:,i));
%     D(:,i) = abs(D(:,i));
%     D(:,i) = D(:,i) - min(D(:,i));
    
    
    % Energy Threshold
    De = nip_energy(D(:,i));
    idx = find(De < 0.1*max(De));
    De(idx) = 0;
    for j = 1:3
        idx3 = (idx-1)*3+j;
        D(idx3,i) = 0;
    end    
    D(:,i) = D(:,i)./norm(D(:,i));
 
    
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
    
    finalD = [finalD,finalDi(:,2:end)];
end

for i=1:size(finalD,2)
   finalD(:,i) = finalD(:,i)./norm(finalD(:,i)); 
end

% Eliminate Redundant elements of the dictionary
% Df = [];
% while true
%     if ~isempty(Df)
%         Dfe = nip_energy(Df);
%     else
%         Dfe = [];
%     end
%     finalDe = nip_energy(finalD);
%     finalDe = aff*finalDe;
%     
%     Rho = corr([Dfe,finalDe]);
%     
%     rho = Rho(:,size(Dfe,2)+1);
%     rho(size(Dfe,2)+1) = 0;
%     idx = find(rho > 0.1);
%     if ~isempty(idx)
%         Df = [Df,finalD(:,1)+sum(finalD(:,idx),2)];
%         finalD(:,[1,idx]) =[];
%     else
%         Df = [Df,finalD(:,1)];
%         finalD(:,1) = [];
%     end
%     
%     if isempty(finalD)
%         break
%     end
% end
% for i = 1:size(Df,2)
%     Df(:,i) = Df(:,i)./norm(Df(:,i));
% end
% L = model.L*Df;
L = model.L*finalD;

if verbose
    figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
    Nh = ceil(sqrt(Np));
    Nw = ceil(sqrt(Np));
    ha = tight_subplot(Nh, Nw,  0.05, 0.1, 0.1);
end

for i = 1:size(finalD,2)
    nip_reconstruction3d(model.cortex,finalD(:,i),opt3d)
end

%% Compute hyperparameters and show results
h = nip_kalman_hyper(model.y,L);
if verbose
    figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
    plot(model.t,h')
    title('Temporal evolution of the hyperparameters')
    pause(0.01)
end