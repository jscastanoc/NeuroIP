% Get results from real data in the cluster

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 27 Aug 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization and creation of the structures used in some functions
close all; clc; clear

addpath('/home/jscastanoc/bbci/toolbox/startup')
addpath('/home/jscastanoc/svn_test/matgrid')
addpath('/home/jscastanoc/NeuroIP')

addpath('../../external/source_toolbox/haufe/')
addpath('../../external/source_toolbox/nolte/')
addpath('../../external/source_toolbox/simulations/')
rng('default')
rng('shuffle');
nip_init();
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
% cfg.fs = 100; % Sample frequency (Hz)
% cfg.t = 0:1/cfg.fs:1.5; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa model;


startup_bbcicluster;
% session_list = get_session_list('D:\bbci\investigation\projects\AudioVisualSpeller');
session_list = get_session_list('/home/jscastanoc/NeuroIP/scripts_papers/SFLEX-TFMxNE/');

stim_target = [31:36 81:86 111:116];
stim_nontarget = [11:16 61:66 91:96];


% Conditions = 3;
Conditions = 2;

jobs_c = 1;
% methods = {'LOR','TF-MxNE','S+T','S-FLEX'};
methods = {'S+T'};
for c_meth = 1:numel(methods)
    
    copy_res = {};
    cur_jobs = [];
    for icond = 1:Conditions
%     for icond = 2
        for ii = 1:length(session_list)
%         for ii = 11
                fname = strcat(session_list{ii}, '/AudiVisual_Depend_', num2str(icond) , '*');
                [ epo, epo_r] = ...
                    stdERPanalysis(fname, {stim_target stim_nontarget; 'Target' 'Non-Target'}, ...
                    'hp_filt', [.4, .2, 3, 30], 'lp_filt', [17 25 3 50],'varReject_freq_band', [5 25],  'ref_ival', []);
                epo = proc_selectChannels(epo, scalpChannels);
                epo_avg = proc_average(epo);
                dmy = proc_selectClasses(epo_avg, 'Non-Target');
                jobs(jobs_c) =  mgsub({'J_rec', 'te'},'testbench_realdata_cluster_norm', ...
                    {methods{c_meth}, dmy}, 'qsub_opts', '-l h_vmem=8G');
                cur_jobs = [cur_jobs jobs(jobs_c)];
                jobs_c = jobs_c + 1;
                dir = strcat('/home/jscastanoc/Results_final_real_norm/',session_list{ii});
                file_name = strcat(dir,'/',methods{c_meth},'Cond',num2str(icond),'.mat');
                copy_res{end+1} = file_name;
        end        
    end
    mgwait(cur_jobs);
    % Copy results
    aux = 1;
    for cc_res = cur_jobs
        or = strcat('/home/jscastanoc/svn_test/matgrid/jobs/',num2str(cc_res),'/mgjob_results.mat');
        dest = copy_res{aux};
        try
            movefile(or,dest);
        catch
            err_var = true;
            save(file_name,'err_var');
        end
        aux = aux+1;
    end
    mgclear(cur_jobs);
end