% Script to generate dataset in the cluster and then rename it and copy it
% in to another directory.

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 19 Aug 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization and creation of the structures used in some functions
close all; clc; clear;
addpath('/home/jscastanoc/svn_test/matgrid');
addpath('/home/jscastanoc/bbci/toolbox/startup');
addpath('/home/jscastanoc/NeuroIP');
nip_init();
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
cfg.fs = 120; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:1.5; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa;
Nd = size(model.cortex.vc,1);

Ntrials = [5 20 50 100 250];
act_sources = [1 3 5];
% Ntrials = [210];
snr_meas = 0;
snr_bio = 5;
% snr_bio = [5]
snr = [];
% t = model.t;
% L = model.L;
% fs = model.fs
% cortex = model.cortex;
% clear model;

n_exp = 1;

distmat = graphrbf(model.cortex);

jobs_c = 1;
methods = {'LOR','TF-MxNE','S+T','S-FLEX'};
% methods = {'S-FLEX'};
dir_base = '/home/jscastanoc/simulated/montreal_sampleall_false/';
dir_results = '/home/jscastanoc/results/montreal_sampleall_false/';
dir_error = '/mnt/data/error/montreal_sampleall_false/';

for c_meth = 1:numel(methods)
    for j = 1:n_exp
        cur_jobs = [];
        copy_res = {};
        for l = 1:length(act_sources)
            for i = Ntrials
                dir = strcat(dir_base,num2str(act_sources(l)),'/');
                file_name = strcat(dir,'Exp',num2str(j),'Ntrials',...
                    num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                load(file_name);
                model.y = y;
                
                jobs(jobs_c) = mgsub({'J_rec', 'time', 'er'},'solvers_ip', ...
                    {model , methods{c_meth}, Jclean}, 'qsub_opts', '-l h_vmem=8G');
                cur_jobs = [cur_jobs jobs(jobs_c)];
                jobs_c = jobs_c + 1;
                
                dir = strcat(dir_results,num2str(act_sources(l)));
                file_name = strcat(dir,'/',methods{c_meth},'Exp',num2str(j),'Ntrials',...
                    num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                copy_res{end+1} = file_name;
                
                dir = strcat(dir_error,num2str(act_sources(l)));
                file_name = strcat(dir,'/',method,'Exp',num2str(j),'Ntrials',...
                    num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                error_file{end+1}= file_name;
                
            end
        end
        mgwait(cur_jobs);
        
        % Copy results
        aux = 1;
        for cc_res = cur_jobs
            or = strcat('/home/jscastanoc/svn_test/matgrid/jobs/',num2str(cc_res),'/mgjob_results.mat');
            try
                load(or);
                J_rec = mgjob.results{1};
                time = mgjob.results{2};
                
                dest = copy_res{aux};
                save(dest,'J_rec','time','dir_base');
                
                dest = error_file{aux};
                save(dest,'er')
                delete(or);
            catch
                err_var = true;
                save(dest,'err_var');
            end
            aux = aux+1;
        end
        mgclear(cur_jobs);
    end
end