function [J_rec,te] = testbench_realdata_cluster(methods,dmy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Testbench for LORETA SFLEX TFMxNE and Proposed %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 26 Aug 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%
nip_init();
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
    basis{1} = basis{1}/sum(basis{1}(:));

end

% Depth compensation
Lsloreta = nip_translf(model.L);
Winv = full(sloreta_invweights(Lsloreta));
Winv = cat(3,Winv(1:end/3,1:end/3),Winv(end/3 +1:2*end/3,end/3 +1:2*end/3),Winv(2*end/3 +1 :end,2*end/3 +1:end));
for i = 1:3
    Lsloreta(:,:,i) = Lsloreta(:,:,i)*Winv(:,:,i);
end
model.L = nip_translf(Lsloreta);

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
        
    case 'S+T'
        options.iter = 50;
        options.spatial_reg =1.5;
        options.temp_reg = 0.5;
        options.tol = 1e-2;
        [J_rec,~] = nip_sflex_tfmxne(model.y,model.L,basis{1},options);
    otherwise
        error(strcat('Nah! ',methods{i},' is not available'))
end
J_est = nip_trans_solution(J_est);
for i = 1:3
    J_est(:,:,i) = Winv(:,:,i)*J_est(:,:,i);
end
J_est = nip_trans_solution(J_est);
te = toc;