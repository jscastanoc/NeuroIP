% Script for non uniform priors.
% Juan S. Castano C. 4 Mar 2013
clear; close all; clc;
nip_init();
Nd = 4000;
load(strcat('data/montreal',num2str(Nd),'_10-10.mat'))
% load(strcat('data/default_data'))
% cortex_mesh = vol;

cfg.L = L;
cfg.cortex = cortex_mesh;
cfg.fs = 200;
cfg.t = 0:1/cfg.fs:1;

model = nip_create_model(cfg);
clear cfg;

[Laplacian QG] = nip_neighbor_mat(model.cortex);
Nt = length(model.t);

A    = triangulation2adjacency(model.cortex.faces);
D   = compute_distance_graph(A);

%%
close all; clc;
sigma = 2;
aff = exp(-D/sigma);
% plot(aff(1,:),'r');
% hold on;
% plot(aff(:,1),'b');



act = exp(-(model.t-0.75).^2/0.05).*sin(2*pi*10*model.t);
x = nip_simulate_activity(model.cortex,Laplacian, [15 10 30], ...
        act,model.t);
[~,idx] =  max(sum(x.^2));

x = aff*x;

figure('Units','normalized','Position',[0.1 0.1 0.3 0.3])
nip_reconstruction3d(model.cortex,aff(:,idx).*sum(x.^2,2),gca);