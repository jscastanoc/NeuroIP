% Plot script for Figure 1 Real Data

close all, clc, clear;

addpath('../../external/source_toolbox/haufe/')
addpath('../../external/source_toolbox/nolte/')
addpath('../../external/source_toolbox/simulations/')
rng('default')
rng('shuffle');
nip_init();

global DATA_DIR;
DATA_DIR = '/mnt/data/Datasets/';
startup_bbci;
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

session_list = get_session_list(strcat(fileparts(which('script_figure_1'))));

cond = 1;
% for j = 9:15 %9 13 14 7
% for j = [7 9 13 14]
stim_target = [31:36 81:86 111:116];
stim_nontarget = [11:16 61:66 91:96];





for j = 9
    session_n = session_list{j};
    for i = 1:numel(methods)
        close all
        ff = figure('Units','normalized','position',[0.1 0.1 0.33 0.6]);
        fname = strcat(session_list{j}, '/AudiVisual_Depend_', num2str(cond) , '*');
        [ epo, epo_r, ivals, cnt, misc ] = stdERPanalysis(fname, {stim_target stim_nontarget; 'Target' 'Non-Target'}, ...
            'hp_filt', [.4, .2, 3, 30], 'lp_filt', [17 25 3 50],'varReject_freq_band', [3 25],  'ref_ival', []);
        
        res_dir = '/mnt/data/ResultadosBerlin/Results_final_real/';
        file_name = strcat(res_dir,session_n,'/',methods{i},'Cond',num2str(cond),'.mat');
        load(file_name);
        J_rec = mgjob.results{1};
        epo = proc_selectChannels(epo, scalpChannels);
        epo_avg = proc_average(epo);
        % data of non-target
        dmy = proc_selectClasses(epo_avg, 'Non-Target');
        
        t0or = 210;
        t1or = 450;
        
        t0 = t0or-dmy.t(1);
        t0 = t0/1000;
        t0 = round(t0*dmy.fs);
        t1 = t1or-dmy.t(1);
        t1 = t1/1000;
        t1 = round(t1*dmy.fs);
        subplot(3,3,1:3)
        plot(dmy.t,dmy.x,'k');
        xlim([dmy.t(1) dmy.t(end)]);
        ylim([min(dmy.x(:)) max(dmy.x(:))]);
        hold on
        plot([t0or t0or],[min(dmy.x(:)), max(dmy.x(:))],'r')
        plot([t1or t1or],[min(dmy.x(:)), max(dmy.x(:))],'r')
        xlabel('Time (ms)')
        ylabel('Amp. (uV)')
        % Options for 3d plot
        options3d.thres = 0;
        options3d.colormap = 'jet';
        options3d.crange = [0 0.18];
        
        subplotxl(3,3,4)
        mnt = getElectrodePositions(dmy.clab);
        scalpPlot(mnt, dmy.x(t0,:), 'scalePos', 'none','colAx',[min(dmy.x(:)),max(dmy.x(:))]);
        subplotxl(3,3,5)
        options3d.view= [90 0];
        nip_reconstruction3d(model.cortex, J_rec(:,t0),options3d);
        subplotxl(3,3,6)
        options3d.view= [-90 0];
        nip_reconstruction3d(model.cortex, J_rec(:,t0),options3d);
        
        subplotxl(3,3,7)
        mnt = getElectrodePositions(dmy.clab);
        scalpPlot(mnt, dmy.x(t1,:), 'scalePos', 'none','colAx',[min(dmy.x(:)),max(dmy.x(:))]);
        subplotxl(3,3,8)
        options3d.view= [90 0];
        nip_reconstruction3d(model.cortex, J_rec(:,t1),options3d);
        subplotxl(3,3,9)
        options3d.view= [-90 0];
        nip_reconstruction3d(model.cortex, J_rec(:,t1),options3d);
    end
    fig_name =strcat('/mnt/data/ResultadosBerlin/Figures1/Session',num2str(j),'Cond',num2str(cond));
    % savefig(fig_name, ff, 'pdf');
    savefig(fig_name, ff, 'eps');
    % savefig(fig_name, ff, 'png');
end