close all, clc, clear;

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

% temp = sa.V_cortex10K;
temp = sa.V_cortex_coarse;
L = nip_translf(temp); % Leadfield matrix
L = L(find(ismember(clab_example,clab)),:);
% L = nip_translf(temp); % Leadfield matrix
clear temp

% cfg.cortex = sa.cortex10K;
cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = 100; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:1.5; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa cfg;

act_s = 1;
exper = 10;
Ntrials = 100;
BioNoise = -5;

sim_dir = 'D:/Datasets/sim_trials_final/';
file_name = strcat(sim_dir,num2str(act_s),'/Exp',...
    num2str(exper),'Ntrials',num2str(Ntrials),...
    'BioNoise',num2str(BioNoise),'.mat');
load(file_name);


model.y = mgjob.results{1};
J_clean = mgjob.results{2};
actidx = mgjob.results{3};
% clear mgjob;
methods = {'LOR','S-FLEX','TF-MxNE','S+T'};

for i = 1:numel(methods)
    res_dir = 'D:/Results_final/';
    file_name = strcat(res_dir,num2str(act_s),'/',methods{i},'Exp',...
        num2str(exper),'Ntrials',num2str(Ntrials),...
        'BioNoise',num2str(BioNoise),'.mat');
    load(file_name);
    J_rec = mgjob.results{1};
    ff(i) = figure('Units','normalized','position',[0.1 0.1 0.3 0.7]);
    subplot(4,1,1:2)
    nip_reconstruction3d(model.cortex, sqrt(sum(J_rec.^2,2)),gca);
    colorbar off    
    subplot(4,1,3)
    plot(model.t, J_rec((actidx-1)*3 +1,:));
    legend('1','2','3')
    subplot(4,1,4)
    real_act = J_clean((actidx-1)*3 +1,:);
    real_act = real_act./repmat(max(abs(real_act),[],2),1,model.Nt);
    rec_act = J_rec((actidx-1)*3 +1,:);
    rec_act = rec_act./repmat(max(abs(rec_act),[],2),1,model.Nt);
    ERR = sqrt((real_act - rec_act).^2);
    imagesc(ERR)
    colormap jet
end

% ff(end + 1) = figure('Units','normalized','position',[0.1 0.1 0.3 0.7]);
% J_rec = J_clean;
% J_rec(1,:)=1e-1*ones(1,model.Nt);
% subplot(4,1,1:2)
% nip_reconstruction3d(model.cortex, sqrt(sum(J_rec.^2,2)),gca);
% hold on
% scatter3(model.cortex.vc(actidx,1),model.cortex.vc(actidx,2),model.cortex.vc(actidx,3),'filled');
% colorbar off
% subplot(4,1,3)
% plot(model.t, J_rec((actidx-1)*3 +1,:));
% legend('1','2','3')

