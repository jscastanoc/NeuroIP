% Script to make some test on correlated sources
% Juan S. Castano
% 29 jan 2013;

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
% Laplacian = eye(model.Nd);


act = [sin(2*pi*10*model.t) ; sin(2*pi*10*model.t) ]; 
x = nip_simulate_activity(model.cortex,Laplacian, QG, [4 -3 7; 4 3 7], ...
        act,model.t);
model.y = model.L*x;


model.cortex.vertices = model.cortex.vertices -repmat(mean(model.cortex.vertices),model.Nd,1)
model.cortex.vertices(:,2) = -model.cortex.vertices(:,2)
a =10*(2*pi)/360;
rotx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
model.cortex.vertices = 10*(rotx*model.cortex.vertices')';

figure(1)
subplot(2,4,1)
% nip_reconstruction3d(model.cortex,mean(x.^2,2),gca);
nip_glassbrain(model.cortex.vertices, mean(x.^2,2));

dip_act = [];
L = model.L;


% Front

regions{1} = ~cellfun('isempty',regexpi(elec.label,'(fp|af|f)([13579]|z)'));
regions{2} = ~cellfun('isempty',regexpi(elec.label,'(fp|af|f)([2468]|z|10)'));
regions{3} = ~cellfun('isempty',regexpi(elec.label,'(ft|fc|c|t|tp|cp)([2468]|z|10)'));
regions{4} = ~cellfun('isempty',regexpi(elec.label,'(ft|fc|c|t|tp|cp)([13579]|z)'));
regions{5} = ~cellfun('isempty',regexpi(elec.label,'(p|po|o|i)([2468]|z|10)'));
regions{6} = ~cellfun('isempty',regexpi(elec.label,'(p|po|o|i)([13579]|z)'));

tmp = regions;
clear regions;
regions{1} = tmp{1}|tmp{3}|tmp{5};
regions{2} = tmp{2}|tmp{4}|tmp{6};

for j = 1:numel(regions)
Q{j} = sparse([1,model.Nd]);
ch = regions{j};
L = model.L(ch,:);
for i = 1:model.Nd
    delta(i) = 1/(L(:,i)'*L(:,i));
end
InvCov = spm_inv(model.y(ch,:)*model.y(ch,:)');            
for i = 1:model.Nd
        Q{j}(i) = 1/(L(:,i)'*InvCov*L(:,i));
        Q{j}(i) = Q{j}(i)/delta(i);
end

Q{j} = Q{j}/max(abs(Q{j}));
subplot(2,4,j+1)
nip_glassbrain(model.cortex.vertices, Q{j});
end

% x_rec = nip_loreta(model.y,model.L,100*diag(Q{1}+Q{2}));
% subplot(2,5,j+2)
% nip_glassbrain(model.cortex.vertices, mean(x_rec.^2,2))

figure
for i = 1:numel(regions)
    subplot(3,2,i)
    scatter3(elec.chanpos(regions{i},1),elec.chanpos(regions{i},2),elec.chanpos(regions{i},3),'xr')
    hold on
    scatter3(cortex_mesh.vertices(:,1),cortex_mesh.vertices(:,2),cortex_mesh.vertices(:,3))
    
    axis equal
end
% 
% % selected_ch = [3 5 7 10 11 12 13 20 21 22 24 26 29 30];
% % s_ch = ismember((1:model.Nc),selected_ch);
% % [~,selected_ch] = find(~s_ch);
% selected_ch = 1:model.Nc
% 
% seg = 1:model.Nd;
% aux = model.L;
% model.L = model.L(selected_ch,seg);
% model.cortex.vertices = model.cortex.vertices(seg,:);
% Q = Q(seg,seg);
% model.Nd = length(seg);
% 
% model.y = model.y(selected_ch,:)
% x_rec = nip_loreta(model.y, model.L, Q);
% % x = x(seg,:);
% 
% subplot(3,2,2);
% % nip_reconstruction3d(model.cortex,mean(x_rec.^2,2),gca);
% nip_glassbrain(model.cortex, mean(x_rec.^2,2));
% % nip_reconstruction3d(model.cortex,diag(Q),gca);
% 
% subplot(3,2,3);
% imagesc(x*x');
% % stem(diag(x*x'))
% colorbar;
% model.L = aux;
% sum(diag(model.L*x*x'*model.L'))
% 
% subplot(3,2,4);
% imagesc(x_rec*x_rec');
% sum(diag( model.L(selected_ch,seg)*x_rec*x_rec'* model.L(selected_ch,seg)'))
% colorbar;
% 
% subplot(3,2,5)
% plot(model.y')
% 
% subplot(3,2,6)
% plot((model.L*x_rec)')