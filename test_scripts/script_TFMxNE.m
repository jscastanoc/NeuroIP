% function main()
% Script Tf-MxNE
clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%
% Simulation %
%%%%%%%%%%%%%%

% Number of dipoles(to use with the provided data. It should not be
% difficult to adapt it to your own data
Nd = 1000; % 1000, 2000, 4000, 8000

% load data
load(strcat('../data/nocons_or/montreal',num2str(Nd),'_full1shell.mat'))

cfg.L = lf; % Leadfield matrix
cfg.cortex = cortex_mesh; % Cortical mesh
cfg.fs = 200; % Sample frequency
cfg.t = 0:1/cfg.fs:0.20; % Time vector (seconds)
cfg.elec = elec; %Struct with labels (field: label) and position (field: elecpos and chanpos) of the channels

% Put everything together in the same struct
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec



% Activity that's going to be simulated each row corresponds to the activity in one dipole.
act(1,:) = -model.t.^2.*exp(-model.t/0.02).*sin(2*pi*10*model.t); % Actividad a simular
act(1,:) = act(1,:)/max(abs(act(1,:)));

% For example, if you want other 2 active dipoles
x = 10*(model.t-3); % Some time scaling, ignore this
x = model.t*40-5;
act(2,:) = exp(-(x).^2).*(-x.^5+x.^2);
act(3,:) = -exp(-(model.t-0.15).^2/0.005).*sin(2*pi*5*model.t-3.5);


% J-> Simulated activity on the brain
[J, idx] = nip_simulate_activity(model.cortex.vertices,[15 20 15; 5 0 15; 15 -20 15 ], ...
        act,randn(size(act,1),3),model.t);
% [J, idx] = nip_simulate_activity(model.cortex, [5 0 15], ...
%         act,model.t);

% gaussian blobs in the cortical surface:
fuzzy = nip_fuzzy_sources(model.cortex,0.01);
index = (1:3:model.Nd);

% At this point theres only completly focal sources, we can apply a low
% pass filter by using the gaussian blobs obtained above
for i = 0:2
    J(index+i,:) = fuzzy*J(index+i,:); 
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is how the simulation looks like %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),gca);
rec_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
subplot(1,2,1)
plot(J(:,:)')
hold on



% Simulate EEG with the simulated brain activity.
model.y = model.L*J;
model.y = nip_addnoise(model.y,10);

Lcomp = nip_depthcomp(model.L,0.3);

[J_rec,~] = nip_tfmxne(model.y,Lcomp,...
    struct('spatial_reg',0.5,'temp_reg',0.05));

subplot(1,2,2)
plot(J_rec(:,:)');
figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);
title('TFMxNE')


