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

Ntrials = [5 20 50 100 250];
% Ntrials = [100];
% act_sources = [1 3 5];
Nact = [1 3 5];
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

Nexp = [1:25];

jobs_c = 1;
% methods = {'LOR','TF-MxNE','STOUT','S-FLEX','KAL','IRA3','IRA5','LOR_PROJ'};
methods = {'S-FLEX','STOUT','TF-MxNE'};
dir_base = '/home/jscastanoc/simulated/montreal_sampleall_false/';
dir_results = '/home/jscastanoc/results/montreal_sampleall_false/';
dir_error = '/home/jscastanoc/error/montreal_sampleall_false/';
% dir_base = '/home/jscastanoc/simulated/montreal_sampleall_false/';
% dir_results = '/home/jscastanoc/results/montreal_sampleall_false/depth/';
% dir_error = '/home/jscastanoc/error/montreal_sampleall_false/depth/';
% dir_base = '/mnt/data/Master_Results/Datasets/simulated/montreal_sampleall_false/';
% dir_results = '/mnt/data/results/montreal_sampleall_false/';
% dir_error = '/mnt/data/error/montreal_sampleall_false/';
depth = {'none','sLORETA','Lnorm'};
depth ={'Lnorm'};
tic
for l = 1:numel(depth)
    cur_jobs = [];
    copy_res = {};
    error_file = {};
    for i = Nexp
        for c_meth = 1:numel(methods)
            for j = Nact
                for k = Ntrials
                    method = methods{c_meth};
                    dir = strcat(dir_base,num2str(j),'/');
                    file_name = strcat(dir,'Exp',num2str(i),'Ntrials',...
                        num2str(k),'BioNoise',num2str(snr_bio),'.mat');
                    load(file_name);
                    load_data;
                    model.fs = fs;
                    model.y = nip_addnoise(y,8);
                    model.Nt = size(y,2);
                    model.t = 0:1/fs:model.Nt/fs;
                    
                    resnorm = norm(model.y - model.L*Jclean, 'fro')/norm(model.y, 'fro')
                    
                    
                    jobs(jobs_c) = mgsub({'J_rec', 'time', 'er', 'extra'},'thesis_core', ...
                        {Nexp,Nact,Ntrials,method,dir_base,dir_results,...
                        dir_error,snr_bio,model, Jclean, actidx,depth{l},resnorm}, 'qsub_opts', '-l h_vmem=4G');
                    cur_jobs = [cur_jobs jobs(jobs_c)];
                    jobs_c = jobs_c + 1;
                    
                    dir = strcat(dir_results,num2str(j));
                    file_name = strcat(dir,'/',methods{c_meth},'Exp',num2str(i),'Ntrials',...
                        num2str(k),'BioNoise',num2str(snr_bio),depth{l},'.mat');
                    copy_res{end+1} = file_name;
                    
                    dir = strcat(dir_error,num2str(j));
                    file_name = strcat(dir,'/',method,'Exp',num2str(i),'Ntrials',...
                        num2str(k),'BioNoise',num2str(snr_bio),depth{l},'.mat');
                    error_file{end+1}= file_name;
                    
                end
            end
        end
        
    end
    
    mgwait(cur_jobs);
    
    % Copy results
    aux = 1;
    for cc_res = cur_jobs
        or = strcat('/home/jscastanoc/svn_test/matgrid/jobs/',num2str(cc_res),'/mgjob_results.mat');
        load(or);
        J_rec = mgjob.results{1};
        time = mgjob.results{2};
        er = mgjob.results{3};
        extra = mgjob.results{4};
        
        dest = copy_res{aux};
        save(dest,'J_rec','time','dir_base','extra');
        
        dest = error_file{aux};
        save(dest,'er')
        delete(or);
        aux = aux+1;
    end
    mgclear(cur_jobs);
end
toc