% Plot script for Figure 2 Real Data

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

methods = {'S+T'};

startup_bbcicluster;
session_list = get_session_list('D:\bbci\investigation\projects\AudioVisualSpeller');

cond = 2;
% for j = 9:15 %9 13 14 7
% for j = [7 9 13 14]
stim_target = [31:36 81:86 111:116];
stim_nontarget = [11:16 61:66 91:96];

    

J_avg = zeros(size(model.L,2),1);
class = 'Non-Target';
n_total = 1;
for j = 1:15
    session_n = session_list{j};
    for i = 1:numel(methods)
        res_dir = 'D:/Results_final_real/';
        file_name = strcat(res_dir,session_n,'/',methods{i},'Cond',num2str(cond),'.mat');
        load(file_name);
        J_rec = mgjob.results{1};
        
    end
    J_rec = J_rec/max(abs(J_rec(:)));
    J_rec = sqrt(sum(J_rec.^2,2));
    if ~sum(isnan(J_rec))
        J_avg = J_avg + J_rec;
        n_total = n_total + 1;
    end
end
J_avg = J_avg/n_total;

% Options for 3d plot
options3d.thres = 0;
options3d.colormap = 'jet';
% options3d.crange = [0 1];

ff = figure('Units','normalized','position',[0.1 0.1 0.3 0.2]);
subplotxl(1,2,1)
options3d.view= [90 0];
nip_reconstruction3d(model.cortex, J_avg,options3d);
subplotxl(1,2,2)
options3d.view= [-90 0];
nip_reconstruction3d(model.cortex, J_avg,options3d);




fig_name =strcat('D:/Figures2/', class,'Cond',num2str(cond));
savefig(fig_name, ff, 'pdf');
savefig(fig_name, ff, 'eps');
savefig(fig_name, ff, 'png');