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
labels = clab_example{find(ismember(clab_example,clab))};

% cfg.cortex = sa.cortex10K;
cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = 100; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:1.5; % Time vector (seconds)
model = nip_create_model(cfg);
clear L cfg;

methods = {'LOR','TF-MxNE','S+T'};
methods = {'S+T'};

startup_bbcicluster;
session_list = get_session_list('D:\bbci\investigation\projects\AudioVisualSpeller');

cond = 2;
% for j = 9:15 %9 13 14 7
% for j = [7 9 13 14]
stim_target = [31:36 81:86 111:116];
stim_nontarget = [11:16 61:66 91:96];

for j = 1
    session_n = session_list{j};
for i = 1:numel(methods)
    fname = strcat(session_list{j}, '/AudiVisual_Depend_', num2str(cond) , '*');

% [epo, epo_r] = stdERPanalysis(fname, {stim_target stim_nontarget; 'Target' 'Non-Target'}, ...
%     'plotting', 1, 'hp_filt', [.4, .2, 3, 30], 'lp_filt', [17 25 3 50], ...
%     'disp_ival', [-200 800],'varReject_freq_band', [5 25],  'ref_ival', []);
% break
[ epo, epo_r, ivals, cnt, misc ] = stdERPanalysis(fname, {stim_target stim_nontarget; 'Target' 'Non-Target'}, ...
    'hp_filt', [.4, .2, 3, 30], 'lp_filt', [17 25 3 50],'varReject_freq_band', [5 25],  'ref_ival', []);
    res_dir = 'D:/Results_final_real/';
    file_name = strcat(res_dir,session_n,'/',methods{i},'Cond',num2str(cond),'.mat');
    load(file_name);
    J_rec = mgjob.results{1};
    ff(i) = figure('Units','normalized','position',[0.1 0.1 0.3 0.4]);
%     subplot(4,1,1:2)
%     nip_reconstruction3d(model.cortex, sqrt(sum(J_rec.^2,2)),gca);
%     nip_reconstruction3d(model.cortex, J_rec(:,16),gca);
    colorbar off    
%     subplot(4,1,3)
%     plot(model.t, J_rec((actidx-1)*3 +1,:));
%     legend('1','2','3')
%     subplot(4,1,4)
%     real_act = J_clean((actidx-1)*3 +1,:);
%     real_act = real_act./repmat(max(abs(real_act),[],2),1,model.Nt);
%     rec_act = J_rec((actidx-1)*3 +1,:);
%     rec_act = rec_act./repmat(max(abs(rec_act),[],2),1,model.Nt);
%     ERR = sqrt((real_act - rec_act).^2);
%     imagesc(ERR)
%     colormap jet
end
close all
end

