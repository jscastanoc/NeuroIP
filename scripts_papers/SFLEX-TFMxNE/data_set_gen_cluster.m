% Script to generate dataset in the cluster and then rename it and copy it
% in to another directory.

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 19 Aug 2013           %
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
clear L sa;


Ntrials = [5 20 50 100 250];
% Ntrials = [10];
snr_meas = 0;
snr_bio = [-15 -5 5];
% snr_bio = [5]
Nspurious = 200;
snr = [];
t = model.t;
L = model.L;
fs = model.fs
cortex = model.cortex;
clear model;


phase_shift{1} = [0.75];
phase_shift{2} = [0.375 0.75 1.25];
phase_shift{3} = [0.250 0.500 0.750 1.000 1.250];
n_exp = 50;


jobs_c = 1;
for k = 1:length(snr_bio)
    for j = 1:n_exp  
        cur_jobs = [];
        copy_res = {};
        for l = 1:numel(phase_shift)
            ps= phase_shift{l}+ 0.05*randn(1,length(phase_shift{l})); % Phase shift for the sources (in seconds)
            Nact = length(phase_shift);
            fc_wl = 9*ones(1,length(ps)) + 2.5*randn(1,length(ps)); % Central frequency for the wavelet
            source_act = [];
            for m = 1:length(ps)
                source_act = [source_act ;aux_wavelet(fc_wl(m),t,ps(m))];
            end
            snr = [];
            for i = Ntrials
%                 [y, Jclean, actidx] = nip_simtrials(L, cortex.vc, ...
%                     source_act, t, Nspurious , i, snr_meas, snr_bio(k));
                jobs(jobs_c) = ...
                    mgsub({'y', 'Jclean', 'actidx'},'nip_simtrials', ...
                    {L, cortex.vc, source_act, t, Nspurious , i, snr_meas, snr_bio(k)});
                cur_jobs = [cur_jobs jobs(jobs_c)];
                jobs_c = jobs_c + 1;
                dir = strcat('/home/jscastanoc/sim_trials_final/',num2str(length(ps)));
                file_name = strcat(dir,'/Exp',num2str(j),'Ntrials',...
                        num2str(i),'BioNoise',num2str(snr_bio(k)),'.mat');
%                 copy_res{jobs_c} = {file_name,jobs(jobs_c)}
                copy_res{end+1} = file_name;
%                 y_clean = L*Jclean;
%                 noise = y-y_clean;
%                 %             y = model.y;
%                 snr = [snr 20*log10(norm(y_clean)/norm(noise))];
            end
        end        
        mgwait(cur_jobs);
        
        % Copy results
        aux = 1;
        for cc_res = cur_jobs
            or = strcat('/home/jscastanoc/svn_test/matgrid/jobs/',num2str(cc_res),'/mgjob_results.mat');
            dest = copy_res{aux};
            movefile(or,dest);
            aux = aux+1;
        end
        mgclear(cur_jobs);
        
    end
end