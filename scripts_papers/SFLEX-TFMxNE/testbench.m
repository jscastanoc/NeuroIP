%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Testbench for LORETA SFLEX TFMxNE and Proposed %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 13 Aug 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization and creation of the structures used in some functions
close all; clc; clear

addpath('../../external/source_toolbox/haufe/')
addpath('../../external/source_toolbox/nolte/')
addpath('../../external/source_toolbox/simulations/')
rng('default')
rng('shuffle');

warning off

load clab_example;
load clab_10_10;
clab = clab_10_10;
data_name = 'icbm152b_sym';
% data_name = 'montreal';

sa = prepare_sourceanalysis(clab, data_name);

temp = sa.V_cortex10K;
L = nip_translf(temp); % Leadfield matrix
clear temp

cfg.cortex = sa.cortex10K;
cfg.L = L;
cfg.fs = 200; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:2; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa;


%% Generation of the Morlet wavelet %%

phase_shift = [0.5 1.5] ; % Phase shift for the sources (in seconds)
Nact = length(phase_shift);
fc_wl = 5; % Central frequency for the wavelet
for i = 1:Nact
    f0 = fc_wl;
    
    % Normalization terms for the wavelet
    sigma_f = f0/7;
    sigma_t = 1/(2*pi*sigma_f);
    
    % "source" contains the time courses of active sources
    source(i,:) =  real(exp(2*1i*pi*f0*model.t).*...
        exp((-(model.t-phase_shift(i)).^2)/(2*sigma_t^2)));
end
source_act = source;
clear source;

%% Simulation of the EEG (Ntrial Averages)
Ntrials = 50;
snr_meas = 10;
snr_bio = 0;
Nspurious = 200;
[model.y, Jclean, active] = nip_simtrials(model.L, model.cortex.vc, ...
    source_act, model.t, Nspurious , Ntrials, snr_meas, snr_bio);


%% Source reconstruction

methods = {'LOR','S-FLEX','TF-MxNE','S+T'};
s_fig = figure('Units','normalized','position',[0.1 0.1 0.45 0.3]);
t_fig = figure('Units','normalized','position',[0.1 0.1 0.8 0.3]);
e_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);

iter = 1;
for i = 1:numel(methods)
    switch methods{i}
        case 'LOR'
            [Laplacian, ~] = nip_neighbor_mat(model.cortex);
            Q = inv(Laplacian*Laplacian');
            [J_rec,~] = nip_loreta(model.y,model.L,diag(Q));
            if (show_res == true) && iter ==1
                figure(s_fig)
                subplot(numel(methods),1,i)
                nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);
                title(methods{i})
                figure(t_fig)
                subplot(numel(methods),1,i)
                plot(model.t,J_rec(active,:)')
                title(methods{i})
            end
            
        case 'S-FLEX'
            [S, out] = sflex_cortical_dal(model.y, nip_translf(model.L),model.cortex,...
                struct('eps',0.1));
            J_rec = nip_trans_solution(S);
            if (show_res == true) && iter ==1
                figure(s_fig)
                subplot(numel(methods),1,i)
                nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);
                title(methods{i})
                figure(t_fig)
                subplot(numel(methods),1,i)
                plot(model.t,J_rec(active,:)')
                title(methods{i})
            end
            
        case 'TF-MxNE'
            options.iter = 50;
            options.spatial_reg = 4;
            options.temp_reg = 1.2;
            options.tol = 1e-2;
            [J_rec,~] = nip_tfmxne_port(model.y,model.L,options);
            if (show_res == true) && iter ==1
                figure(s_fig)
                subplot(numel(methods),1,i)
                nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);
                title(methods{i})
                figure(t_fig)
                subplot(numel(methods),1,i)
                plot(model.t,J_rec(active,:)')
                title(methods{i})
            end
            
        case 'S+T'
            options.iter = 50;
            options.spatial_reg = 4;
            options.temp_reg = 1.2;
            options.tol = 1e-2;
            [J_rec,~] = nip_sflex_tfmxne(model.y,model.L,basis,options);
            if (show_res == true) && iter ==1
                figure(s_fig)
                subplot(numel(methods),1,i)
                nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);
                title(methods{i})
                figure(t_fig)
                subplot(numel(methods),1,i)
                plot(model.t,J_rec(active,:)')
                title(methods{i})
            end
            
        otherwise
            error(strcat('Nah! ',methods{i},' is not available'))
    end
end

