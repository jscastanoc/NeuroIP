clc, close all, clear;

rng('default')
rng('shuffle');
nip_init();
ltfatstart; % Init Time frequency toolbox
warning off

load clab_example
load clab_10_10;
clab = clab_10_10;
% data_name = 'icbm152b_sym';
data_name = 'montreal';
sa = prepare_sourceanalysis(clab, data_name);

temp = sa.V_cortex_coarse;
L = nip_translf(temp); % Leadfield matrix
L = L(find(ismember(clab_example,clab)),:);
clear temp

cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = 100; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:1.5; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa cfg;


% Generation of the Morlet wavelet (used for the time courses of the active
% dipoles
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

% Simulation of the EEG (Ntrial Averaged)
Ntrials = 50;
snr_meas = 10; % SNR at sensor level
snr_bio = 0; % SNR at source level
Nspurious = 200; % Number of active "spurious" dipoles
[model.y, J_clean, active] = nip_simtrials(model.L, model.cortex.vc, ...
    source_act, model.t, Nspurious , Ntrials, snr_meas, snr_bio);


model.L = nip_depthcomp(model.L,0.2); % Depth bias compensation for the lead field matrix

% Construction of the spatial basis functions
nbasis = size(model.L,2)/3;
iter_basis = [1]; % Width of the spatial basis functions
basis = [];
n = 1;
group = [];
for i = iter_basis
    fuzzy = nip_fuzzy_sources(model.cortex,i);
    basisn = fuzzy(:,randi([1,model.Nd/3],nbasis,1));
    basisn = basisn/norm(basisn(:),1);
    basis{n} = basisn;
    group = [group n*ones(1,nbasis)];
    basis{n} = nip_blobnorm(basis{n},group,struct('norm',1,'norm_group',false)); % Normalization of the spatial basis functions
    n = n+1;
end


options.iter = 50; % Max number of iteration
options.spatial_reg = 0.9;
options.temp_reg = 0.2;
options.tol = 1e-2;
[J_rec,~] = nip_sflex_tfmxne(model.y,model.L,basis{1},options);


%%%%%%%%%%%%%%%%%
% Visualization %
%%%%%%%%%%%%%%%%%
% Simulation - Temporal
figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
plot(model.t,J_clean')
xlabel('Time')
ylabel('Amplitude')

% Simulation - Spatial
figure('Units','normalized','position',[0.2 0.2 0.14 0.14]);
nip_reconstruction3d(model.cortex, sqrt(sum(J_clean.^2,2)), struct('axes',gca));
hold on
scatter3(model.cortex.vc(active,1),model.cortex.vc(active,2),model.cortex.vc(active,3),'ok','filled');


%%% Reconstruction - Temporal
figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
plot(model.t,J_rec')
xlabel('Time')
ylabel('Amplitude')

% Reconstruction - Spatial
figure('Units','normalized','position',[0.2 0.2 0.15 0.2]);
nip_reconstruction3d(model.cortex, sqrt(sum(J_rec.^2,2)), struct('axes',gca));
hold on
h = scatter3(model.cortex.vc(active,1),model.cortex.vc(active,2),model.cortex.vc(active,3),'ok','filled');
