function [] = testbench(methods,iter,show_results,save_results)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Testbench for LORETA SFLEX TFMxNE and Proposed %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 13 Aug 2013           %
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
clab = clab_10_10;
% data_name = 'icbm152b_sym';
data_name = 'montreal';

sa = prepare_sourceanalysis(clab, data_name);

% temp = sa.V_cortex10K;
temp = sa.V_cortex_coarse;
L = nip_translf(temp); % Leadfield matrix
L = L(find(ismember(clab_example,clab)),:);
% L = nip_translf(temp); % Leadfield matrix
clear temp

% cfg.cortex = sa.cortex10K;
cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = 100; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:1.5; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa cfg;


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
    iter_basis = [1];
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

model.L = nip_depthcomp(model.L,0.2);

% show_results = true;
% save_results = false;
use_pre_sim = true;
Ntrials = [100];

if (~save_results && ...
        show_results && ...
        length(Ntrials ==1) && ...
        iter ==1)
    s_fig = figure('Units','normalized','position',[0.1 0.1 0.45 0.3]);
    t_fig = figure('Units','normalized','position',[0.1 0.1 0.8 0.3]);
    e_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
end


for k = Ntrials;
    for j = 1:iter
        if use_pre_sim
            dir = 'D:/Datasets/sim_trials_final/1/';
            file_name = strcat(dir,'Exp',num2str(j),'Ntrials',...
                num2str(k),'BioNoise',num2str(-15),'.mat');
            load(file_name);
            model.y = mgjob.results{1};
            Jclean = mgjob.results{2};
            actidx = mgjob.results{3};
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
                    options.spatial_reg = 0.9;
                    options.temp_reg = 0.2;
                    options.tol = 1e-2;
                    [J_rec,~] = nip_tfmxne_port(model.y,model.L,options);
                    
                case 'S+T'
                    options.iter = 50;
                    options.spatial_reg = 0.9;
                    options.temp_reg = 0.2;
                    options.tol = 1e-2;
                    [J_rec,~] = nip_sflex_tfmxne(model.y,model.L,basis{1},options);
                otherwise
                    error(strcat('Nah! ',methods{i},' is not available'))
            end
            if save_results
                dir = 'D:/Datasets/Results/';
                J_rec = sparse(J_rec);
                file_name = strcat(dir,'Exp',num2str(j),'Ntrials',...
                    num2str(k),'BioNoise',num2str(-5),methods{i},'.mat');
                te = toc;
                save(file_name,'J_rec','te');
            end
            
            if (~save_results && ...
                    show_results && ...
                    length(Ntrials ==1) && ...
                    iter ==1)
                figure(s_fig)
                subplot(1,numel(methods)+1,i+1)
                nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);
                
                figure(t_fig)
                subplot(1,numel(methods)+1,i+1)
                idx = [];
                for i = 1:length(actidx)
                    idx(:,i) = ((actidx(i)-1)*3+1:(actidx(i)-1)*3+3);
                end
                idx = idx(:);
                plot(model.t,J_rec(idx-1,:))
            end
            
        end
        if (~save_results && ...
                show_results && ...
                length(Ntrials ==1) && ...
                iter ==1)
            figure(s_fig)
            subplot(1,numel(methods)+1,1)
            nip_reconstruction3d(model.cortex,sqrt(sum(Jclean.^2,2)),gca);
            
            figure(t_fig)
            subplot(1,numel(methods)+1,1)
            idx = [];
            for i = 1:length(actidx)
                idx(:,i) = ((actidx(i)-1)*3+1:(actidx(i)-1)*3+3);
            end
            idx = idx(:);
            plot(model.t,Jclean(idx-1,:))
        end
        
    end
    
end
