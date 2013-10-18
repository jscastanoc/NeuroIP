% Script to generate dataset using the parallel computing toolbox

%%%%%%%%%%%%%%%%%%%%%%%%%
% Juan S. Castano C.    %
% jscastanoc@gmail.com  %
% 19 Sep 2013           %
%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear
warning off
nip_init();

% Load leadfield data (leadfield, electrode labels etc..)
load_montreal_1010;

cfg.cortex = sa.cortex_coarse;
cfg.L = L;
cfg.fs = 120; % Sample frequency (Hz)
cfg.t = 0:1/cfg.fs:1.5; % Time vector (seconds)
model = nip_create_model(cfg);
clear L sa;


Ntrials = [5 20 50 100 250]; % Number of trials to simulate (results are averaged across trials)
snr_meas = 0; % Signal to noise ratio at the sensor level
snr_bio = [5]; % Signal to noise ratio at source level
Nspurious = 1000; % Number of source with spurious activity in each simulation


% Mean of the time shift of the wavelets (1 3 and 5 active sources)
phase_shift{1} = [0.75];
phase_shift{2} = [0.375 0.75 1.25];
phase_shift{3} = [0.250 0.500 0.750 1.000 1.250];

% Number of Experiments
n_exp = 100;

options.Ntrials = Ntrials;
options.snr_bio = snr_bio;
options.n_exp = n_exp;
options.Nspurious = Nspurious;
options.snr_meas = snr_meas;

sched = findResource('scheduler', 'configuration', 'local');
% sched = parcluster();
job = createJob(sched);
dir_base = '/mnt/data/Master_Results/Datasets/simulated/montreal_sampleall_false/';
for l = 1:numel(phase_shift)
    options.phase_shift = phase_shift{l};
%     parallel_core(model, dir_base, options);
    tasks{l} = createTask(job, @parallel_core, 0, {model, dir_base, options});        
end
submit(job);
% wait(job);
waitForState(job, 'finished');
errmsgs = get(job.Tasks, {'ErrorMessage'});
nonempty = ~cellfun(@isempty, errmsgs);
celldisp(errmsgs(nonempty));
results = getAllOutputArguments(job); 
destroy(job);

copyfile('*.m',dir_base);