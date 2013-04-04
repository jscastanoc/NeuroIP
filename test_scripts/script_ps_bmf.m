% Script to benchmark the iterative regulariation algorithm
% Juan S. Castano C. 26 Jan 2013;
clear; close all; clc;
nip_init();

% Load data (Leadfield, etc...) to be used
Nd = 8000;
load(strcat('data/montreal',num2str(Nd),'_10-10.mat'))
% load(strcat('data/default_data'))
% cortex_mesh = vol;

cfg.L = L;
cfg.cortex = cortex_mesh;
cfg.fs = 200;
cfg.t = 0:1/cfg.fs:1;

model = nip_create_model(cfg);
[Laplacian QG] = nip_neighbor_mat(model.cortex);

Nt = length(model.t);

act1 = exp(-(model.t-0.25).^2/0.05).*sin(2*pi*10*model.t);
act2 = exp(-(model.t-0.75).^2/0.05).*sin(2*pi*10*model.t)
act = [act1; act2];
% figure
% plot(cfg.t,act');
% xlabel('Time (Sec)')
% ylabel('Amplitude (a.u.)')
% pause
% savefig('times_2dip',gcf,'eps')

% figure
% sim = 0.9*ones(size(model.L,2),1);
% sim(1) = 1;
% sim(2) = 0;
% nip_reconstruction3d(model.cortex,sim,gca);
% colormap gray
% colorbar off;
% pause
% savefig('brain_2dip',gcf,'pdf')
% break

[~,index_t1] = max(act(1,:));
[~,index_t2] = max(act(2,:));

x = nip_simulate_activity(model.cortex,Laplacian, [8 -17 30; 15 10 30], ...
        act,model.t);
    
clean_y = model.L*x;


figure
ha = tight_subplot(2,2, 0.01, 0.01, 0.01);
% axes(ha(1))
% nip_reconstruction3d(model.cortex,x(:,index_t1),gca);
% axes(ha(2))
% nip_reconstruction3d(model.cortex,x(:,index_t2),gca);


snr = 5;
model.y = nip_addnoise(clean_y, snr);


Q = inv(Laplacian'*Laplacian);
[x, extras] = nip_loreta(model.y, model.L, Q);

model.J_rec = x;
model.extras_out = extras;

axes(ha(1))
nip_reconstruction3d(model.cortex,x(:,index_t1),gca);
axes(ha(2))
nip_reconstruction3d(model.cortex,x(:,index_t2),gca);


clear Q;
for i = 1:model.Nd
    delta(i) = 1/(L(:,i)'*L(:,i));
end
InvCov = spm_inv(model.y*model.y');            
for i = 1:model.Nd
        Q(i) = 1/(L(:,i)'*InvCov*L(:,i));
        Q(i) = Q(i)/delta(i);
end
Q = diag(Q);
[x, extras] = nip_loreta(model.y, model.L, Q);
axes(ha(3))
nip_reconstruction3d(model.cortex,x(:,index_t1),gca);
axes(ha(4))
nip_reconstruction3d(model.cortex,x(:,index_t2),gca);
pause
savefig('recLOR_BMF',gcf,'pdf')
