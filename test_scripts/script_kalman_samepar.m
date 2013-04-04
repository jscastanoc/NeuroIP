clear; close all; clc;
nip_init();

save_fig = false;
Nd = 2000;
load(strcat('data/montreal',num2str(Nd),'_full.mat'))
% load(strcat('data/default_data'))
% cortex_mesh = vol;

cfg.L = L;
cfg.cortex = cortex_mesh;
cfg.fs = 200;
cfg.t = 0:1/cfg.fs:1;

model = nip_create_model(cfg);
[Laplacian QG] = nip_neighbor_mat(model.cortex);

Nt = length(model.t);
Nc = model.Nc;
Q = zeros([Nd,1]);
Q(600) = 1;
a = [1 -0.5];
x = zeros([Nd,Nt]);
for k = 3:Nt
    for i = 1:Nd
       x(i,k) = a(1)*x(i,k-1) + a(2)*x(i,k-2) + randn*Q(i); 
    end
    if (k >= Nt/2)
%         a = [1.3 -1];
    end
end
x = nip_simulate_activity(cfg.cortex,speye(Nd), [30 30 30], x(600,:), model.t);
figure
nip_reconstruction3d(cfg.cortex, sum(x.^2,2),gca)

y = L*x + 0.0005*randn([Nc, Nt]);
% ,[1 0 -0.5 1 sqrt(0.0005)]
[x_est extra] = nip_kalmanwh(y,L,speye(Nd));

figure
nip_reconstruction3d(cfg.cortex, sum(x_est.^2,2),gca)