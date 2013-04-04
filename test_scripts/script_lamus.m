% Script implementing Solution given by (Lamus et al 2012)
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

F = speye(model.Nd);
Pt{1} = 1*speye(model.Nd);
xt(:,1) = zeros([model.Nd,1]);
alpha = 2;
beta = 1e-10;
for i=1:20
for k = 2: model.Nt
    x_ap = F*xt(:,k-1);
    P_ap = F*Pt{k-1}*F'+diag(Q);
    K = P_ap*model.L'/(model.L*P_ap*model.L'+R);
    xt(:,k) = x_ap+K*(model.y(:,k)-model.L*x_ap);
    Pt{k} = (eye(model.Nd)-K*model.L)*P_ap;
end
for k = model.Nt-1:1
   J(:,k) = P{k}*F'*P{k+1};
   x(:,k) = xt(:,k) + J(:,k)*(x(:,k+1)-xt(:,k+1));
   P{k} = Pt{k}+J(:,k)*(P{k+1}-Pt{k+1})*J(:,k)';
   lagP{k} = P{k}*J(:,k-1);
end
for k=2:model.Nt
   A1 = A1 + (P{k}+x(:,k)*x(:,k)');
   A2 = A2 + (lagP+x(:,k)*x(:,k-1)');
   A3 = A3 + (P{k-1} + x(:,k-1)*x(:,k-1)');
end
A = A1-A2*F'-F*A2'+F*A3*F';
for i = 1:model.Nd
   Q(i) = A(i,i) + 2 *beta;
   Q(i) = Q(i)/(model.Nt+2*(alpha+1));
end
Sigma0 = P{1};
nip_reconstruction3d(model.cortex,mean(x.^2,2),gca);
end


% [x_rec alpha] = nip_loreta(model.y, model.L, diag(Q));
% [x_rec alpha] = nip_loreta(model.y, model.L, inv(Laplacian'*Laplacian));
nip_reconstruction3d(model.cortex,mean(x_rec.^2,2),gca);
% 
% [x_rec alpha] = nip_loreta(model.y, model.L, Laplacian);
% 
% nip_reconstruction3d(model.cortex,diag(Q),gca);