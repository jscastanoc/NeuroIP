function [J_rec,te] = testbench_realdata_cluster(methods,dmy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Testbench for LORETA SFLEX TFMxNE and STOUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 26 Aug 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%
nip_init();
% Initialization and creation of the structures used in some functions
% close all; clc; clear

addpath('external/source_toolbox/haufe/')
addpath('external/source_toolbox/nolte/')
addpath('external/source_toolbox/simulations/')
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
clear temp

% cfg.cortex = sa.cortex10K;
cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = dmy.fs; % Sample frequency (Hz)
cfg.t = dmy.t; % Time vector (seconds)
model = nip_create_model(cfg);
model.y = dmy.x';
model.y = model.y(find(ismember(clab,sa.clab_electrodes)),:);
clear L sa cfg;


depth = 'Lnorm'; % it can be none, Lnorm or sLORETA-based depth compensation
switch depth
    case 'none'
        L = model.L;
        Winv = [];
    case 'Lnorm'
        gamma = 0.7; % How strong the depth compensation is?
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;
end
clear extras;
model.L = L;

transM = (1+(1/model.Nc))*eye(model.Nc)-(1/model.Nc)*ones(model.Nc,model.Nc);

model.y = transM*model.y;
model.L = transM*model.L;


if sum(ismember(methods,{'STOUT','S-FLEX'}))
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
    basis{1} = basis{1}/sum(basis{1}(:));
    
end


switch methods
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
        options.spatial_reg = 1.5;
        options.temp_reg = 0.5;
        options.tol = 1e-2;
        [J_rec,~] = nip_tfmxne_port(model.y,model.L,options);
        
    case 'STOUT'
        sigma = 1; % Width of the gaussian bells or cortical blobs used as spatial dictionary
        
        % This function can be used to create spatial dictionaries using a forward
        % model taking into account only the cortex surface (model.cortex has the
        % 3d graph representing the cortical surface)
        B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
        
        % Normalize the spatial basis functions
        B = nip_blobnorm(B,'norm',2);
        
                
        % Options for the inversion
        % The ratio between spatial_reg and temp_reg depends on the snr of the EEG
        % However, in general a ratio of 1:3 should work ok.
        spatial_reg = 80 ; % Sparsity in the spatial domain
        temp_reg =  10; % Sparsity in the time-frequency domain
        gof = 0.85;
        a = 10;  %  Time shift for the Short Time Fourier Transform (STFT).
        m = 100; %Frequency bins for the STFT.
        lipschitz = [];
        % By setting 'optimgof' to true, the regularization parameters will be
        % modified to get the desired gof. If false, then the user-select reg.
        % parameters are used for the solution
        [J_est, extras] = nip_stout(model.y, L, B,'optimgof',true,...
            'sreg',spatial_reg,'treg',temp_reg,'gof', gof, 'a',a ,'m',m,...
            'lipschitz', lipschitz,'Winv',Winv);
    otherwise
        error(strcat('Nah! ',methods{i},' is not available'))
end
te = toc;