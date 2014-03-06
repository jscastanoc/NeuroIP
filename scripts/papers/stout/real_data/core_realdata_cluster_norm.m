function [J_rec,te] = testbench_realdata_cluster(method,dmy)

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
clab = dmy.clab;
data_name = 'montreal';

sa = prepare_sourceanalysis(clab, data_name);

% Reduce lead field to the available sensors
temp = sa.V_cortex_coarse;
L = nip_translf(temp); 
L = L(find(ismember(clab_example,sa.clab_electrodes)),:);
clear temp


cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = dmy.fs; % Sample frequency (Hz)
cfg.t = dmy.t; % Time vector (seconds)
model = nip_create_model(cfg);
model.y = dmy.x';
model.y = model.y(find(ismember(clab,sa.clab_electrodes)),:);
clear L sa cfg;

depth = 'Lnorm';
switch depth
    case 'none'
        L = model.L;
        Winv =[];
    case 'Lnorm'
        gamma = 0.6;
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;
end
clear extras;

% Reference
transM = eye(model.Nc)-(1/model.Nc)*ones(model.Nc);

model.y = transM*model.y;
L = transM*L;


tic;
resnorm = 0.7;
switch method
    case 'LOR'
        [J_rec, extras] = nip_loreta(model.y,L,'cov',speye(size(L,2)),'Winv',Winv);
    case 'S-FLEX'
        % Spatial dictionary
        sigma = 1;
        B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
        B = nip_blobnorm(B,'norm',2);
        % Options for the inversion
        reg_par = 100;
        [J_rec, extras] = nip_sflex(model.y, L, B,'optimres',true,'regpar', reg_par,'resnorm',resnorm,'Winv',Winv);
    case 'TF-MxNE'
        spatial_reg = 160;
        temp_reg =  10;
        [J_rec, extras] = nip_tfmxne_python(model.y, L, 'optimres',true,'sreg',spatial_reg,'treg',temp_reg,'resnorm',resnorm,'Winv',Winv);
    case 'STOUT'
        % Spatial dictionary
        sigma = 1;
        B = nip_fuzzy_sources(model.cortex, sigma, struct('save',1,'dataset','montreal'));
        B = nip_blobnorm(B,'norm',2);
        
        % Options for the inversion
        spatial_reg = 160;
        temp_reg =  10;
        [J_rec, extras] = nip_stout_python(model.y, L, B,'optimres',true,'sreg',spatial_reg,'treg',temp_reg,'resnorm',resnorm,'Winv',Winv);
    otherwise
        error('%s not available as method',method)
end
te = toc;
