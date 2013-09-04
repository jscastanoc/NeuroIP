%  Compute error meas cluster

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
act_sources = [1 3 5];
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

n_exp = 10;


jobs_c = 1;
methods = {'LOR','TF-MxNE','S+T','S-FLEX'};
errors = [1 2];
dM = graphrbf(model.cortex);
dummy_counter = 0;
total_iter = length(errors)*numel(methods)*(n_exp)*length(act_sources)*length(Ntrials);
for c_er = errors
    for c_meth = 1:numel(methods)
        for j = 1:n_exp
            cur_jobs = [];
            copy_res = {};
            for l = 1:length(act_sources)
                for i = Ntrials
                    dir = strcat('D:/Results_final/',num2str(act_sources(l)));
                    file_name = strcat(dir,'/',methods{c_meth},'Exp',num2str(j),'Ntrials',...
                        num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                    %                 copy_res{jobs_c} = {file_name,jobs(jobs_c)}
                    copy_res{end+1} = file_name;
                    load(file_name)
                    J_rec = mgjob.results{1};
                    dir = strcat('D:/sim_trials_final/',num2str(act_sources(l)));
                    file_name = strcat(dir,'/Exp',num2str(j),'Ntrials',...
                        num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                    load(file_name);
                    J_clean = mgjob.results{2};
                    switch c_er
                        case 1
                            sig1 = sqrt(sum(J_rec.^2,2));
                            sig2 = sqrt(sum(J_clean.^2,2));                            
                            dir = strcat('D:/Errors1/',num2str(act_sources(l)));
                        case 2
                            sig1 = sqrt(sum(J_rec*J_clean',2));
                            sig2 = sqrt(sum(J_clean*J_clean',2));
                            dir = strcat('D:/Errors2/',num2str(act_sources(l)));
                    end
                    sig1 = nip_energy(sig1);
                    sig1 = abs(sig1)/sum(abs(sig1(:)));
                    
                    sig2 = nip_energy(sig2);
                    sig2 = abs(sig2)/sum(abs(sig2(:)));
                    idx = find(sig2);
                    
%                     distance_c = nip_emd(sig2(idx),sig1,dM(:,idx));    
                    distance_c = emd_hat_mex_nes(sig1, sig2(idx),dM(
                    file_name = strcat(dir,'/',methods{c_meth},'Exp',num2str(j),'Ntrials',...
                        num2str(i),'BioNoise',num2str(snr_bio),'.mat');
                    save(file_name,'distance_c');
                    fprintf('Completed %d \n',100*dummy_counter/total_iter);
                    dummy_counter = dummy_counter + 1;
                end
            end
        end
    end
end