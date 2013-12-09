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

load_data;
clear L sa;
Nd = size(model.cortex.vc,1);

% Ntrials = [5 20 50 100 250];
Ntrials = 250;
% act_sources = [1 3 5];
act_sources = [5];
% Ntrials = [210];
snr_meas = 0;
snr_bio = -5;
% snr_bio = [5]
snr = [];
% t = model.t;
% L = model.L;
% fs = model.fs
% cortex = model.cortex;
% clear model;

Nexp = [1:1:1];

jobs_c = 1;
methods = {'LOR','TF-MxNE','STOUT','S-FLEX'};
% methods = {'S-FLEX'};
dir_base = '/home/jscastanoc/simulated/montreal_sampleall_false/';
dir_results = '/home/jscastanoc/results/montreal_sampleall_false/';
dir_error = '/home/jscastanoc/error/montreal_sampleall_false/';
% dir_base = '/mnt/data/Datasets/simulated/montreal_sampleall_false/';
% dir_results = '/mnt/data/results/montreal_sampleall_false/';
% dir_error = '/mnt/data/error/montreal_sampleall_false/';


depth = 'none';
switch depth
    case 'none'
        L = model.L;
    case 'Lnorm'
        gamma = 0.6;
        [L, extras] = nip_depthcomp(model.L,struct('type',depth,'gamma',gamma));
        Winv = extras.Winv;
    case 'sLORETA'
        [L, extras] = nip_depthcomp(model.L,struct('type',depth));
        Winv = extras.Winv;
        
end
clear extras;

for j = 1:n_exp
    cur_jobs = [];
    copy_res = {};
    error_file = {};
    for c_meth = 1:numel(methods)
        for l = 1:length(act_sources)
            for i = Ntrials
                method = methods{c_meth};
                dir = strcat(dir_base,num2str(act_sources(l)),'/');
                file_name = strcat(dir,'Exp',num2str(j),'Ntrials',...
                    num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                load(file_name);
                model.y = y;
                model.gof = gof;
                
                jobs(jobs_c) = mgsub({'J_rec', 'time', 'er'},'solvers_ip', ...
                    {model , methods{c_meth}, Jclean, act_sources(l), actidx}, 'qsub_opts', '-l h_vmem=8G');
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
            
            if ~strcmp(depth,'none')
                J_est = nip_translf(J_rec');
                for i = 1:3
                    J_est(:,:,i) = (Winv(:,:,i)*J_est(:,:,i)')';
                end
                J_rec = nip_translf(J_est)';
            end
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