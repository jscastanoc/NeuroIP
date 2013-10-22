% Simulations to compare sLORETA and norm of L for depth bias compensation
clc; clear; close all;

Nexp = [1:20,32];
Nact = [1];
Ntrials = [250];
methods = {'LOR','STOUT'};
Nparallel = 2;
bionoise = 5;
sim_dir = '/mnt/data/Master_Results/Datasets/simulated/montreal_sampleall_false/';
res_dir = '/mnt/data/Master_Results/STOUT/bin/depthcompensation/sLORETA/';
% res_dir = '/mnt/data/Master_Results/STOUT/bin/depthcompensation/Lnorm/';
% res_dir = '/mnt/data/Master_Results/STOUT/bin/depthcompensation/nonorm/';
err_dir = '/mnt/data/Master_Results/STOUT/error/depthcompensation/sLORETA/';
% err_dir = '/mnt/data/Master_Results/STOUT/error/depthcompensation/Lnorm/';
% err_dir = '/mnt/data/Master_Results/STOUT/error/depthcompensation/nonorm/';
run_parallel_exp(Nexp,Nact,Ntrials,methods,Nparallel, sim_dir,res_dir, err_dir, bionoise);