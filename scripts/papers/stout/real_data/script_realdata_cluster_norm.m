% Get results from real data in the cluster

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 27 Aug 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear

% Path BBCI
addpath('/home/jscastanoc/bbci/toolbox/startup')

% Path matgrid
addpath('/home/jscastanoc/svn_test/matgrid')

% Path NeuroIP toolbox
addpath('/home/jscastanoc/NeuroIP')

% Path Guido Nolte and Stefan Haufe imaging toolbox
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
data_name = 'montreal';

sa = prepare_sourceanalysis(clab, data_name);

temp = sa.V_cortex_coarse;

% Reduce lead field to the available sensors
L = nip_translf(temp); 
L = L(find(ismember(clab_example,clab)),:);
clear temp

cfg.cortex = sa.cortex_coarse;
cfg.L = L;
model = nip_create_model(cfg);
clear L sa model;


startup_bbcicluster;

session_list = get_session_list(pwd);

stim_target = [31:36 81:86 111:116];
stim_nontarget = [11:16 61:66 91:96];

% Reconstruct for conditions from 1 ... Conditions
Conditions = 2; 

jobs_c = 1;

% Imaging method and stimulus type ('Target' or 'Non-Target')
methods = {'STOUT'};
type = 'Non-Target';
for c_meth = 1:numel(methods)    
    copy_res = {};
    cur_jobs = [];
    for icond = 1:Conditions
        for ii = 1:length(session_list)
            fname = strcat(session_list{ii}, '/AudiVisual_Depend_', num2str(icond) , '*');
            [ epo, epo_r] = ...
                stdERPanalysis(fname, {stim_target stim_nontarget; 'Target' 'Non-Target'}, ...
                'hp_filt', [.4, .2, 3, 30], 'lp_filt', [17 25 3 50],'varReject_freq_band', [5 25],  'ref_ival', []);
            epo = proc_selectChannels(epo, scalpChannels);
            epo_avg = proc_average(epo);
            dmy = proc_selectClasses(epo_avg, type);
            jobs(jobs_c) =  mgsub({'J_rec', 'te'},'core_realdata_cluster_norm', ...
                {methods{c_meth}, dmy}, 'qsub_opts', '-l h_vmem=4G');
            cur_jobs = [cur_jobs jobs(jobs_c)];
            jobs_c = jobs_c + 1;
            dir = strcat('/home/jscastanoc/real_data_results/',session_list{ii});
            file_name = strcat(dir,'/',methods{c_meth},'Cond',num2str(icond),type,'.mat');
            copy_res{end+1} = file_name;
        end
    end
    mgwait(cur_jobs);
    
    aux = 1;
    
    % Copy results to the final directory
    for cc_res = cur_jobs
        or = strcat('/home/jscastanoc/svn_test/matgrid/jobs/',num2str(cc_res),'/mgjob_results.mat');
        dest = copy_res{aux};
        % movefile(or,dest);
        load(or);
        J_rec = mgjob.results{1};
        time = mgjob.results{2};
        save(dest,'J_rec','time');
        delete(or);
        aux = aux+1;
    end
    mgclear(cur_jobs);
end
