% Script implementing Champagne (Wipf et al 2010)
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


act = [sin(2*pi*10*model.t);sin(2*pi*10*model.t + pi/4)]; 
x = nip_simulate_activity(model.cortex,Laplacian, QG, [-40 0 30; 40 0 30], ...
        act,model.t);
    
model.y = model.L*x;
model.y = nip_addnoise(model.y,-10);



figure
nip_reconstruction3d(model.cortex,mean(x.^2,2),gca);

% Q = model.L'*diag(diag(model.y*model.y'))*model.L;
% Q = Q/(max(diag(Q)));
% Q = inv(Laplacian'*Laplacian);
% Q = model.L'*model.y*model.y'*model.L;
Q = 1+0.001*randn(model.Nd,1);

R = 0.01*eye(model.Nc);



figure 
for j=1:20
for i = 1:model.Nd
   S_b = R + model.L(:,i)*diag(Q(i))*model.L(:,i)';   
   S_bi = inv(S_b);
   X = (Q(i)*model.L(:,i)'*S_bi)*model.y;
   Z(i) = model.L(:,i)'*S_bi*model.L(:,i);   
   Zsqr(i) = sqrt(Z(i));
   Zsqri(i) = 1/Zsqr(i);
   Q(i) = Zsqri(i)*sqrt((Zsqr(i)*X*X'*Zsqr(i)))*Zsqri(i);
end

% nip_reconstruction3d(model.cortex,Q,gca);
%     S_bb = R + model.L*diag(Q)*model.L'; 
%     G_funp = G_fun;
%     G_fun = trace((1/model.Nt)*model.y*model.y'/S_bb)+log(det(S_bb))
%     if abs(G_fun-G_funp) <= 0.1
%         break
%     end
end
Q = Q/max(Q);
[x_rec alpha] = nip_loreta(model.y, model.L, diag(Q));
% [x_rec alpha] = nip_loreta(model.y, model.L, inv(Laplacian'*Laplacian));
nip_reconstruction3d(model.cortex,mean(x_rec.^2,2),gca);
% 
% [x_rec alpha] = nip_loreta(model.y, model.L, Laplacian);
% 
% nip_reconstruction3d(model.cortex,diag(Q),gca);