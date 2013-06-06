function Q = BMTF(y,L)
% Q = BMTF(y,L)
%
% Juan S. Castano
% 25 May 2013

% Time Frequency representation and Parafac decomposition
Factors = TF_PARAFAC(y);
Nfac = size(Factors{1},2);
Nd = size(L,2);

% Reconstruction with each factor
Q = zeros(Nd,Nfac);
for i = 1:Nfac
    Q(:,i)  = nip_lcmv(Factors{1}(:,i),L.^2);
%     figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
%     nip_reconstruction3d(model.cortex,J_rec{i},gca);
end
%
% figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
% nip_reconstruction3d(model.cortex,nip_lcmv(model.y,model.L)',gca);
