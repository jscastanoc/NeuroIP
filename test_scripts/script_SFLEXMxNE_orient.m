function main()
% Script PARAFAC
clear; close all; clc;
nip_init();

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preproceso / simulacion %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Numero de dipolos a considerar
Nd = 2000; % 1000, 2000, 4000, 8000

% Cargar datos(lead field, mesh del cerebro etc...
load(strcat('../data/const_or/montreal',num2str(Nd),'_full.mat'))

leadfield = L;
cfg.L = L;
cfg.cortex = cortex_mesh;
cfg.fs = 300; % Frecuencia de muestreo para la simulacion
cfg.t = 0:1/cfg.fs:0.20; % Vector de tiempo
cfg.elec = elec;

% Crea una estructura (model) con los datos cargados y declarados arriba
model = nip_create_model(cfg);
clear cfg L cortex_mesh eeg_std head elec lf

act(1,:) = -model.t.^2.*exp(-model.t/0.02).*sin(2*pi*10*model.t); % Actividad a simular
act(1,:) = act(1,:)/max(abs(act(1,:)));
x = 10*(model.t-3);
x = model.t*40-5;
% act(2,:) = exp(-(x).^2).*(-x.^5+x.^2);
% act(3,:) = -exp(-(model.t-0.15).^2/0.005).*sin(2*pi*5*model.t-3.5);
% act(1,:) = sin(2*pi*20*model.t+1.5);
% act(2,:) = sin(2*pi*40*model.t);

for i = 1:3
    acti{i} = act; % Actividad a simular
end
act = acti;
clear acti;

% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, 1*[15 20 15; 15 5 5; 15 -20 15; 15 5 25 ], ...
%         act,model.t);
[J, idx] = nip_simulate_activity(model.cortex, [5 0 15], ...
        act,model.t);
nn = (1:3:model.Nd*3);
n = 1;
for i = nn
   JJ(n,:) = J(i,:);
   n = n+1;
end
J = JJ;
clear JJ
% [J, idx] = nip_simulate_activity(model.cortex,Laplacian, size(act,1), ...
%         act,model.t);
fuzzy = nip_fuzzy_sources(model.cortex,0.1);
J = fuzzy*J; % J simulado FINAL


figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
nip_reconstruction3d(model.cortex,sqrt(sum(J.^2,2)),gca);

rec_fig = figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
subplot(1,2,1)
plot(J')
hold on

model.y = model.L*J;
model.y = nip_addnoise(model.y,10);

nbasis = 256;
iter_basis = [0.05 0.1 0.2];
basis = [];
for i = iter_basis
    fuzzy = nip_fuzzy_sources(model.cortex,i);
    basisn = fuzzy(:,randi([1,model.Nd],nbasis,1));
    basis = [basis basisn];
end
basis = reshape(repmat(basis(:)',1,1),model.Nd,nbasis*length(iter_basis));
model.L = 10*model.L*basis;
model.Nd = size(model.L,2);

a = 10; M = 200; L=a*M; 
c = dgtreal(model.y','gauss',a,M);
% c = permute(c,[3 1 2]);
or = idgtreal(c,'gauss',a,M);

T = size(c,2);
K = size(c,1);
Z = sparse(model.Nd,K*T);
Y = Z;


J_rec = sparse(model.Nd,model.Nt);
J_recf = sparse(model.Nd,M);
tau = 1;
mu = 0.01;
lambda = 0;
tempGY = model.L'*model.y;
for i = 1:model.Nd
    lambda = max(norm(tempGY(i,:)),lambda);
end
lambda = lambda*0.002;
clear tempGY;

for i = 1:5
   i
   Z_0 = Z;
   
   [act_dip,~] = ind2sub(size(Y),find(Y));
   if ~isempty(act_dip)
       for j = unique(act_dip)'
          aux = idgtreal(reshape(full(Y(j,:)'),K,T),'gauss',a,M);
          J_rec(j,:) = aux(1:model.Nt);
       end
   end
   error = model.y - model.L*J_rec;
   err_trans = conj(dgtreal(error','gauss',a,M));
   err_trans = permute(err_trans,[3 1 2]);
   
   err_trans = reshape(err_trans,model.Nc,[]);   
   arg_prox = (Y - mu*model.L'*err_trans);
   clear err_trans
   Z = prox(arg_prox,mu,lambda);
   clear arg_prox
   tau_0 = tau;
   tau = (1+sqrt(1+4*tau^2))/2;
   Y = Z + ((tau_0-1)/tau)*(Z-Z_0);
end
[act_dip,~] = ind2sub(size(Z),find(Z));
for j = unique(act_dip)'
    J_recf(j,:) = idgtreal(reshape(full(Z(j,:)),K,T),'gauss',a,M);
end

J_recf = basis*J_recf;

subplot(1,2,2)
plot(J_recf(:,1:model.Nt)');
figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
imagesc(full(abs(Z)))
figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
nip_reconstruction3d(model.cortex,sqrt(sum(J_recf(:,1:model.Nt).^2,2)),gca);
title('TFMxNE')


[J_rec,~] = nip_loreta(model.y,model.L,speye(model.Nd));
J_rec = basis*J_rec;
model.L = leadfield;
[model.Nc model.Nd] = size(model.L);
figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);
title('LORETA')
figure
plot(J_rec')
end

function Z = prox(Y,mu,lambda)

[P,K] = size(Y);

zpk = sparse(P,K);
A = (Y./abs(Y));
A(find(isnan(A))) = 0;
B = max(abs(Y)-lambda,zpk);
AB = A.*B;
clear A B;

parfor p = 1:P
    C(p) = sqrt(sum(max( (abs(Y(p,:))-lambda).^2,0),2));
    aux = mu./C(p);
    aux(find(isnan(aux))) = 0;
    CC(p) = max(1 - aux,0);
end

Z = AB.*repmat(CC',[1,K]);

end
