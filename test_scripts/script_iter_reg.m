% Script to benchmark the iterative regulariation algorithm
% Juan S. Castano C. 26 Jan 2013;
clear; close all; clc;
nip_init();

% Load data (Leadfield, etc...) to be used
Nd = 2000;
load(strcat('data/montreal',num2str(Nd),'_10-10'))

cfg.L = L;
cfg.cortex = cortex_mesh;
cfg.t = 0:1/256:0.5;

model = nip_create_model(cfg);
[Laplacian QG] = nip_neighbor_mat(model.cortex);


act = [sin(2*pi*10*model.t) ; sin(2*pi*10*model.t) ]; 
x = nip_simulate_activity(model.cortex,Laplacian, QG, [-100 0 0;100 0 0], ...
        act,model.t);
    
model.y = model.L*x;


figure
nip_reconstruction3d(model.cortex,mean(x.^2,2),gca);

Q = model.L'*diag(diag(model.y*model.y'))*model.L;
Q = Q/(max(diag(Q)));
Q = inv(Laplacian'*Laplacian);
Q = model.L'*model.y*model.y'*model.L;

tic
[x_rec extras] = nip_iterreg(model.y,model.L,Q, Laplacian, 3);
toc

% 
% [x_rec alpha] = nip_loreta(model.y, model.L, Laplacian);
% 
figure 
nip_reconstruction3d(model.cortex,mean(x_rec.^2,2),gca);
% nip_reconstruction3d(model.cortex,diag(Q),gca);