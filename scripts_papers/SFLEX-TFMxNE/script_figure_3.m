% Plot script for Figure 3 Real Data

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


% for j = 9:15 %9 13 14 7
% for j = [7 9 13 14]
stim_target = [31:36 81:86 111:116];
stim_nontarget = [11:16 61:66 91:96];

    
class1 = 'Non-Target';
n_total = 1;
cond1 = 1;
for j = 1:15
    session_n = session_list{j};
    for i = 1:numel(methods)
        res_dir = 'D:/Results_final_real/';
        file_name = strcat(res_dir,session_n,'/',methods{i},'Cond',num2str(cond1),'.mat');
        load(file_name);
        J_rec = mgjob.results{1};
        
    end
    J_rec = J_rec/max(abs(J_rec(:)));
    J_rec = sqrt(sum(J_rec.^2,2));
    if ~sum(isnan(J_rec))
        J_v(:,n_total) = J_rec;
        n_total = n_total + 1;
    end
end

class2 = 'Non-Target';
n_total = 1;
cond2 = 2;
for j = 1:15
    session_n = session_list{j};
    for i = 1:numel(methods)
        res_dir = 'D:/Results_final_real/';
        file_name = strcat(res_dir,session_n,'/',methods{i},'Cond',num2str(cond2),'.mat');
        load(file_name);
        J_rec = mgjob.results{1};        
    end
    J_rec = J_rec/max(abs(J_rec(:)));
    J_rec = sqrt(sum(J_rec.^2,2));
    if ~sum(isnan(J_rec))
        J_a(:,n_total) = J_rec;
        n_total = n_total + 1;
    end
end

for i = 1:12
J_v1(:,i) = nip_energy(J_v(:,i+2));
J_a1(:,i) = nip_energy(J_a(:,i));
end
[h, pval, ci, stats] = ttest(J_v1',J_a1');
J_avg = stats.tstat;
% J_avg(find(J_avg>0) = J_avg(find(J_avg>0)))
% Options for 3d plot
options3d.thres = 0;
options3d.colormap = 'jet';
% options3d.crange = [0 1];

viewdirs = {[0 -1e-6 -1], [-1 0 0], [0 1 0], [0 1e-6 1], [1 0 0], [0 -1 0]};

S = nip_trans_solution(J_avg);
S = J_avg;
% S = J_avg;
figure; showsurface(model.cortex, struct('myviewdir', viewdirs{1}, 'colorbars', 1, 'myfontsize', 30), S');
zlab = get(h.cb, 'ylabel');
set(h.cb, 'fontsize', 24)
set(zlab, 'String', 'A.U.', 'fontsize', 30)

ff = figure('Units','normalized','position',[0.1 0.1 0.3 0.2]);
subplotxl(1,2,1)
options3d.view= [90 0];
nip_reconstruction3d(model.cortex, J_avg,options3d);
subplotxl(1,2,2)
options3d.view= [-90 0];
nip_reconstruction3d(model.cortex, J_avg,options3d);



fig_name =strcat('D:/Figures3/',num2str(cond1),num2str(cond2),class1,class2);
savefig(fig_name, ff, 'pdf');
savefig(fig_name, ff, 'eps');
savefig(fig_name, ff, 'png');