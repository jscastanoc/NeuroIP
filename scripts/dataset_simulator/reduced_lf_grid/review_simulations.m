% See simulations.

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 19 Sep 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear all

warning off
nip_init();

% Load leadfield data (leadfield, electrode labels etc..)
load clab_example;
load clab_10_10;
clab = clab_10_10;
% data_name = 'icbm152b_sym';
data_name = 'montreal';


% Prepare leadfield data.
sa = prepare_sourceanalysis(clab, data_name);
temp = sa.V_cortex_coarse;
L = nip_translf(temp); % Leadfield matrix
L = L(find(ismember(clab_example,clab)),:);
clear temp

cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = 120; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:1.5; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa;


Ntrials = [5 20 50 100 250]; % Number of trials to simulate (results are averaged across trials)
snr_meas = 0; % Signal to noise ratio at the sensor level
snr_bio = [5]; % Signal to noise ratio at source level
Nspurious = 200; % Number of source with spurious activity in each simulation


% Mean of the time shift of the wavelets (1 3 and 5 active sources)
phase_shift{1} = [0.75];
phase_shift{2} = [0.375 0.75 1.25];
phase_shift{3} = [0.250 0.500 0.750 1.000 1.250];

% Number of Experiments
n_exp = 50;

options.Ntrials = Ntrials;
options.snr_bio = snr_bio;
options.n_exp = n_exp;
options.Nspurious = Nspurious;
options.snr_meas = snr_meas;

sched = findResource('scheduler', 'configuration', 'local');
% sched = parcluster();
job = createJob(sched);
dir_base = '/mnt/data/Master_Results/Datasets/simulated/montreal_sampleall_false/';
figure('Units','normalized','position',[0.2 0.2 0.14 0.4]);
for l = 1:numel(phase_shift)
    ps = phase_shift{l};
    for k = 1:length(snr_bio)
        for j = 1:n_exp            
            for i = Ntrials
                dir = strcat(dir_base,num2str(length(ps)));
                file_name = strcat(dir,'/Exp',num2str(j),'Ntrials',...
                    num2str(i),'BioNoise',num2str(snr_bio(k)),'.mat');
                load(file_name)
                subplot(3,1,1)
                nip_reconstruction3d(model.cortex,sqrt(sum(Jclean.^2,2)),[]);
                hold on;
                scatter3(model.cortex.vc(actidx,1), ...
                    model.cortex.vc(actidx,2),model.cortex.vc(actidx,3),'ok','filled');
                subplot(3,1,2)                
                plot(model.t,y);
                subplot(3,1,3)   
                cla;
                fullidx = (actidx-1)*3 +1;
                for kk = 1:length(fullidx)
                    plot(model.t,Jclean(fullidx(kk):fullidx(kk) + 2,:));
                    hold on
                end
            end
        end
    end
end