
clear; close all; clc;
nip_init();


% Load data (Leadfield, etc...) to be used
Nd = 2000;
load(strcat('../data/nocons_or/montreal',num2str(Nd),'_full1shell.mat'))
% load(strcat('data/default_data'))
% cortex_mesh = vol;
L=lf;
cfg.L = lf;
cfg.cortex = cortex_mesh;
cfg.fs = 200;
cfg.t = 0:1/cfg.fs:0.5;

model = nip_create_model(cfg);
[Laplacian QG] = nip_neighbor_mat(model.cortex);

Nt = length(model.t);

var = 0.05;
act1 = 1*exp(-(model.t-0.12).^2/var).*sin(2*pi*15*model.t);
act2 = 1*exp(-(model.t-0.33).^2/var).*sin(2*pi*15*model.t);
act3 = 1*exp(-(model.t-0.6).^2/var).*sin(2*pi*10*model.t);
act4 = 1*exp(-(model.t-0.8).^2/var).*sin(2*pi*10*model.t);
% act = [act1;act2;act3;act4];
act = [act1;act4];
% act = [act1;act2; act4];

for i = 1:3
    acti{i} = act; % Actividad a simular
end
act = acti;
clear acti;

[~,index_t1] = max(act{1}(1,:));
[~,index_t2] = max(act{1}(2,:));
% [~,index_t3] = max(act(3,:));
% [~,index_t4] = max(act(4,:));

% x = nip_simulate_activity(model.cortex,Laplacian, [10 -17 30; 10 -7 30; 12 0 30; 15 10 30 ], ...
%         act,model.t);
% x = nip_simulate_activity(model.cortex,Laplacian, [-10 -17 30; 10 -7 30; 15 10 30 ], ...
%         act,model.t);
x = nip_simulate_activity(model.cortex, [15 -17 30; 15 10 30 ], ...
        act,model.t);
    
fuzzy = nip_fuzzy_sources(model.cortex,0.1);
index = (1:3:model.Nd);
for i = 0:2
    x(index+i,:) = fuzzy*x(index+i,:); % J simulado FINAL
end 

clean_y = model.L*x;

snr = 10;
model.y = nip_addnoise(clean_y, snr);

fvec = figure('Units','normalized','Position',[0.1 0.1 0.3 0.3]);
hvec = tight_subplot(1,2,[0.05 0.02],[0.1 0.05],[0.08 0.03]);

[Ut St Vt] = svd(model.y'*model.y);
% [Vt St] = eig(model.y'*model.y);
% (Dt(end)+Dt(end-1))/sum(diag(Dt))
axes(hvec(1))
plot(Vt(:,1:2))


[Us Ss Vs] = svd(model.y);
% [Vs Ss] = eig(model.y*model.y');
% (Ds(end)+Ds(end-1))/sum(diag(Ds))
axes(hvec(2))
plot(Vs(:,1:2))


i = 0;
cumvar = 0;
while i < 2
    i = i+1;
%     for j = 1:model.Nc
%         y_ap{i}(j,:) = Vs(j,end-i+1)*Vt(:,end-i+1)';
%     end
    y_ap{i} = Us(:,i)*Ss(i,i)*Vs(:,i)';
    cumvar = cumvar + Ss(i,i)/sum(diag(Ss))
%     cumvar = cumvar + Ss(end-i+1,end-i+1)/sum(diag(Ss))
end
Nk = i;

fmeas = figure('Units','normalized','Position',[0.1 0.1 0.3 0.3]);
hmeas = tight_subplot(1,Nk+1,[0.05 0.02],[0.1 0.05],[0.08 0.03]);
axes(hmeas(1))
nip_sens_meas(model.y,gca);
LQpL = {};
for i = 1:Nk
    axes(hmeas(i+1))
    nip_sens_meas(y_ap{i},gca);
    Qdyn{i} = nip_lcmv(y_ap{i}, model.L);
    Qdyn{i} = Qdyn{i}/max(Qdyn{i});
    LQpL{end+1} = L*diag(Qdyn{i})*L';
end
Qbmf{1} = nip_lcmv(model.y,model.L);

fig_sim = figure('Units','normalized','Position',[0.1 0.1 0.3 0.3]);
% hberry = tight_subplot(2,3,[0.05 0.02],[0.1 0.05],[0.08 0.03]);

nip_reconstruction3d(model.cortex,mean(x.^2,2),gca);

% for i =1:Nk
%     axes(hberry(i+1))    
% %     plot(Q{i})
%     nip_reconstruction3d(model.cortex,Q{i}',gca);
% %     axes(hberry(i+1))
% %     colorbar
% end
    
% w{1} = exp(-(model.t-0.2).^2/var);
% w{2} = exp(-(model.t-0.4).^2/var);
% w{3} = exp(-(model.t-0.6).^2/var);
% w{4} = exp(-(model.t-0.8).^2/var);
% x_rec = zeros(model.Nd, model.Nt);
Qe{1} = eye(model.Nc)/model.Nc;
% for i = 1:Nk
%     Qt(:,i) = Qdyn{i}; 
% end
h = zeros(Nk,model.Nt);
model.y = model.y/max(max(abs(model.y)));
off = 10;
for k = 1:model.Nt-off
%     MVB   = spm_mvb(model.y(:,k),model.L,[],Qt,Qe,16);
    Qk = zeros(model.Nd,1);
    win = k:k+off;
    YY = model.y(:,win)*model.y(:,win)';
    [Cy,h,Ph,F]   = spm_reml_sc(YY,[],[Qe LQpL],1);
%     h(:,k) = diag(MVB.cp);
%     h(:,k) = MVB.h;
    
    hyp(:,k) = h(2:end);
    for i = 1:Nk
        Qk = Qk + hyp(i,k)*Qdyn{i}';
    end
    Qk = Qk/max(Qk);
    [x_dyn(:,k), ~] = nip_loreta(model.y(:,k) , model.L, diag(Qk));

end
[x_bmf, ~] = nip_loreta(model.y , model.L, diag(Qbmf{1}));

figure
plot(hyp');

x_rec = x_dyn;

axes(hberry(2))
cla
nip_reconstruction3d(model.cortex, x_rec(:,index_t1), gca);
axes(hberry(3))
cla
nip_reconstruction3d(model.cortex, x_rec(:,index_t2), gca);
axes(hberry(4))
cla
nip_reconstruction3d(model.cortex, x_rec(:,index_t3), gca);
axes(hberry(5))
cla
nip_reconstruction3d(model.cortex, x_rec(:,index_t4), gca);

times_3d(model.cortex, x_dyn);
times_3d(model.cortex, x_bmf);

% figure
% nip_reconstruction3d(model.cortex, mean(x_dyn.^2,2), gca);
% figure
% nip_reconstruction3d(model.cortex, mean(x_bmf.^2,2), gca);

figure
nip_reconstruction3d(model.cortex, Qdyn{1}', gca);
figure
nip_reconstruction3d(model.cortex, Qdyn{2}', gca);