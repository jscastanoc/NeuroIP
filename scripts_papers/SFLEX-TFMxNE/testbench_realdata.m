function [J_rec,te, model] = testbench_realdata(methods,file_name,show_results,save_results,dmy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Testbench for LORETA SFLEX TFMxNE and Proposed %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 26 Aug 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%
nip_init();
if nargin == 1
    iter = 1;
end

% Initialization and creation of the structures used in some functions
% close all; clc; clear

addpath('../../external/source_toolbox/haufe/')
addpath('../../external/source_toolbox/nolte/')
addpath('../../external/source_toolbox/simulations/')
rng('default')
rng('shuffle');

warning off

load clab_example;
load clab_10_10;
% clab = clab_10_10;
clab = dmy.clab;
% data_name = 'icbm152b_sym';
data_name = 'montreal';

sa = prepare_sourceanalysis(clab, data_name);

% temp = sa.V_cortex10K;
temp = sa.V_cortex_coarse;
L = nip_translf(temp); % Leadfield matrix
L = L(find(ismember(clab_example,sa.clab_electrodes)),:);
% L = nip_translf(temp); % Leadfield matrix
clear temp

% cfg.cortex = sa.cortex10K;
cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = dmy.fs; % Sample frequency (Hz)
cfg.t = dmy.t; % Time vector (seconds)
model = nip_create_model(cfg);
model.y = dmy.x';
model.y = model.y(find(ismember(clab,sa.clab_electrodes)),:);
% model.y = (find(ismember(clab_example,clab)),:);
clear L sa cfg;

transM = (1+(1/model.Nc))*eye(model.Nc)-(1/model.Nc)*ones(model.Nc,model.Nc);

model.y = transM*model.y;
model.L = transM*model.L;


%% Generation of the Morlet wavelet and Simulation of the EEG (Ntrial Averages) %%

% phase_shift = [0.5 1.5] ; % Phase shift for the sources (in seconds)
% Nact = length(phase_shift);
% fc_wl = 5; % Central frequency for the wavelet
% for i = 1:Nact
%     f0 = fc_wl;
%     
%     % Normalization terms for the wavelet
%     sigma_f = f0/7;
%     sigma_t = 1/(2*pi*sigma_f);
%     
%     % "source" contains the time courses of active sources
%     source(i,:) =  real(exp(2*1i*pi*f0*model.t).*...
%         exp((-(model.t-phase_shift(i)).^2)/(2*sigma_t^2)));
% end
% source_act = source;
% clear source;

% Ntrials = 50;
% snr_meas = 10;
% snr_bio = 0;
% Nspurious = 200;
% [model.y, Jclean, active] = nip_simtrials(model.L, model.cortex.vc, ...
%     source_act, model.t, Nspurious , Ntrials, snr_meas, snr_bio);


%% Source reconstruction

% methods = {'LOR','S-FLEX','TF-MxNE','S+T'};





if sum(ismember(methods,{'S+T','S-FLEX'}))
    nbasis = size(model.L,2)/3;
    iter_basis = [1.5];
    basis = [];
    n = 1;
    group = [];
    for i = iter_basis
        fuzzy = nip_fuzzy_sources(model.cortex,i);
        basisn = fuzzy(:,randi([1,model.Nd/3],nbasis,1));
        basisn = basisn/norm(basisn(:),1);
        basis{n} = basisn;
        group = [group n*ones(1,nbasis)];
        n = n+1;
    end
    basis{1} = nip_blobnorm(basis{1},group,struct('norm',1,'norm_group',true));
end

model.L = nip_depthcomp(model.L,0.1);

% show_results = true;
% save_results = false;
use_pre_sim = true;

if (~save_results && ...
        show_results )
    s_fig = figure('Units','normalized','position',[0.1 0.1 0.45 0.3]);
    t_fig = figure('Units','normalized','position',[0.1 0.1 0.8 0.3]);
    e_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
end



for i = 1:numel(methods)
    tic
    switch methods{i}
        case 'LOR'
            %                     [Laplacian] = nip_neighbor_mat(model.cortex);
            %                     Q = inv(Laplacian*Laplacian');
            Q = speye(model.Nd);
            [J_rec,~] = nip_loreta(model.y,model.L,Q);
            
        case 'S-FLEX'
            [S, out] = sflex_cortical_dal(model.y, nip_translf(model.L),basis,...
                struct('eps',0.1));
            J_rec = nip_trans_solution(S);
            
        case 'TF-MxNE'
            options.iter = 50;
            options.spatial_reg = 2.5;
            options.temp_reg = 0.5;
            options.tol = 1e-2;
            [J_rec,~] = nip_tfmxne_port(model.y,model.L,options);
            
        case 'S+T'
            options.iter = 50;
            options.spatial_reg =2.5;
            options.temp_reg = 0.5;
            options.tol = 1e-2;
            [J_rec,~] = nip_sflex_tfmxne(model.y,model.L,basis{1},options);
        otherwise
            error(strcat('Nah! ',methods{i},' is not available'))
    end
    if save_results
        dir = 'D:/Datasets/Results_realdata/res_';
        J_rec = sparse(J_rec);
        file_name = strcat(dir,file_name,'.mat');
        te = toc;
        save(file_name,'J_rec','te');
    end
    
    if (~save_results && ...
            show_results)
        figure(s_fig)
        subplot(1,numel(methods)+1,i+1)
        nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);
        
        figure(t_fig)
        subplot(1,numel(methods)+1,i+1)
        plot(model.t,J_rec)
    end
    te = toc;
end