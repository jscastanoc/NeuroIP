function main()
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
cfg.fs = 300; % Sample frequency
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

% Map it to each of the "directions" i.e. acti{1}->x acti{2}->y acti{3}->z
clear acti;

% J-> Simulated activity on the brain
[J, idx] = nip_simulate_activity(model.cortex,[15 20 15; 5 0 15; 15 -20 15 ], ...
        act,randn(size(act,1),3),model.t);
% [J, idx] = nip_simulate_activity(model.cortex, [5 0 15], ...
%         act,model.t);

% gaussian blobs in the cortical surface:
fuzzy = nip_fuzzy_sources(model.cortex,0.05);
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



% Initialization of the TF-MxNE algorithm
a = 10; M = model.fs; L=a*M; 
c = dgtreal(model.y','gauss',a,M);
T = size(c,2);
K = size(c,1);
Z = sparse(model.Nd,K*T);
Y = Z;
J_rec = sparse(model.Nd,model.Nt);
J_recf = sparse(model.Nd,M);


tau = 1;
% lambda = 0.001; % Time regularization parameter
mu = 0; % Spatial regularization parameter
tempGY = model.L'*model.y;
for i = 1:model.Nd
    basepar = max(norm(tempGY(i,:),2),mu);
end
mu = basepar*0.5; % Spatial regularization parameter
lambda = basepar*0.01; % Time regularization parameter
clear tempGY;

%Normalization of the Lead field matrix to Compensate source depth.
index = (1:3:model.Nd);
gamma = 0.3; % How strong is the depth compensation (0<gamma<1)
for i = index
    norm_term = sqrt((norm(model.L(:,i),2)^2+ ...
        norm(model.L(:,i+1),2)^2 + norm(model.L(:,i+2),2)^2)^gamma);
    model.L(:,i:i+2) = model.L(:,i:i+2)/norm_term;
end

for i = 1:5
    tic
   fprintf('Iteration # %d ',i)
   % line 5 of algorithm 1 (see Ref paper)
   Z_0 = Z; 
   
   
   % line 6 of algorithm 1 (see Ref paper)
   [act_dip,~] = ind2sub(size(Y),find(Y));
   if ~isempty(act_dip)
       for j = unique(act_dip)'
          aux = idgtreal(reshape(full(Y(j,:)'),K,T),'gauss',a,M);
          J_rec(j,:) = aux(1:model.Nt);
       end
   end
   error = model.y - model.L*J_rec;
   err_trans = dgtreal(error.','gauss',a,M);
   err_trans = permute(err_trans,[3 1 2]);  
   err_trans = reshape(err_trans,model.Nc,[]);   
   arg_prox = (Y + mu*model.L'*err_trans);
   clear err_trans
   Z = prox(arg_prox,mu,lambda);
   clear arg_prox
   
   % line 7 of algorithm 1 (see Ref paper)
   tau_0 = tau;
   
   % line 8 of algorithm 1 (see Ref paper)
   tau = (1+sqrt(1+4*tau^2))/2;
   
   % line 9 of algorithm 1 (see Ref paper)
   Y = Z + ((tau_0-1)/tau)*(Z-Z_0);
   
   fprintf('Elapsed time: %f seconds\n',toc)
end
[act_dip,~] = ind2sub(size(Z),find(Z));

fprintf('Transforming solution to the time domain: \n%d non-zero time series \n '...
    , length(unique(act_dip)))
for j = unique(act_dip)'
    J_recf(j,:) = idgtreal(reshape(full(Z(j,:)),K,T),'gauss',a,M);
end
subplot(1,2,2)
plot(J_recf(:,1:model.Nt)');
figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
imagesc(full(abs(Z)))
figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
nip_reconstruction3d(model.cortex,sqrt(sum(J_recf(:,1:model.Nt).^2,2)),gca);
title('TFMxNE')

% [J_rec,~] = nip_loreta(model.y,model.L,speye(model.Nd));
% figure('Units','normalized','position',[0.1 0.1 0.3 0.3]);
% nip_reconstruction3d(model.cortex,sqrt(sum(J_rec.^2,2)),gca);
% title('LORETA')
% figure
% plot(J_rec')

beep
end

function Z = prox(Y,mu,lambda)

[P,K] = size(Y);

zpk = sparse(P,K);

% index = 1:3:P;
% n = 1;
% for p = index
%         normY(n,:) = sqrt(sum(Y(p:p+2,:).^2));
%         n = n+1;
% end
% normY = reshape(repmat(normY(:),3,1),P,K);
% A = Y./normY;

A = Y./abs(Y);
A(find(isnan(A))) = 0;
B = max(abs(Y)-lambda,zpk);
AB = A.*B;
clear A B;


% parfor p = 1:P
%     dip_num = ceil(p/3)
%     C(p) = sqrt(sum(max((sqrt(norm(Y(dip_num,:))^2 + ...
%         norm(Y(dip_num+1,:))^2 + norm(Y(dip_num+2,:))^2)-lambda).^2,0),2));
%     aux = mu./C(p);
%     aux(find(isnan(aux))) = 0;
%     CC(p) = max(1 - aux,0);
% end

parfor p = 1:P
    C(p) = sqrt(sum(max( (abs(Y(p,:))-lambda).^2,0),2));
    aux = mu./C(p);
    aux(find(isnan(aux))) = 0;
    CC(p) = max(1 - aux,0);
end

Z = AB.*repmat(CC',[1,K]);

end

