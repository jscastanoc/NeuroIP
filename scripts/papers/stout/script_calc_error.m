% Script to compute error STOUT

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 24 Sep 2013           %
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


% Ntrials = [5 20 50 100 250]; % Number of trials to simulate (results are averaged across trials)
Ntrials = [5 50 100];
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
options.cortex = model.cortex;
options.L = model.L;
clear model;


sched = findResource('scheduler', 'configuration', 'local');
% sched = parcluster();
job = createJob(sched);
dir_sim = '/mnt/data/Datasets/simulated/montreal_sampleall_false/';
dir_results = '/mnt/data/results/montreal_sampleall_false/';
dir_error = '/mnt/data/error/montreal_sampleall_false/';
method = {'LOR','TF-MxNE','S+T','S-FLEX'};
for l = 1:numel(method)
    tasks{l} = createTask(job, @parallel_error, 0, {method{l}, dir_sim, dir_results, dir_error, options});  
%     parallel_error(method{l}, dir_sim, dir_results, dir_error, options);
end
submit(job);
% wait(job);
waitForState(job, 'finished');
errmsgs = get(job.Tasks, {'ErrorMessage'});
nonempty = ~cellfun(@isempty, errmsgs);
celldisp(errmsgs(nonempty));
results = getAllOutputArguments(job); 
destroy(job);