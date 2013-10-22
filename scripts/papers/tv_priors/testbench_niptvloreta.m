% nip_tvloreta testbench

%% Init
clear; close all; clc;
nip_init();
verbose = true;

addpath('external');
addpath('external/kkmeans');
% Load data, define sample rate, number of eeg channels, etc...
load_data;

% Simulate brain activity and generate pseudo EEG
gen_eeg;

% Show simulated activity
if verbose
    nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),struct('view',[90 0]));
    pause(0.01)
end


% Spatial dictionary in case we solve with, for example, S-FLEX
sp_dict = nip_fuzzy_sources(model.cortex,2,struct('dataset','montreal','save',true));
idx = find(sp_dict < 0.05*max(abs(sp_dict(:))));
sp_dict(idx) = 0;
sp_dict = sparse(sp_dict);

J_est = nip_tvloreta(model.y,model.t,model.L,sp_dict,model.cortex);