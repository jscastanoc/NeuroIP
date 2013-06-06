function Q = ARDTF(y,L, basis)
% Q = SFLEXTF(y,L, basis)
%
% Juan S. Castano
% 25 May 2013

% Time Frequency representation and Parafac decomposition
Factors = TF_PARAFAC(y);
Nfac = size(Factors{1},2);

% SFLEX Preprocessing

nbasis = size(basis,2);
Nd = size(L,2);
Nc = size(L,1);

% Reconstruction with each factor
Q = zeros(Nd,Nfac);
Qe = eye(Nc);
for i = 1:Nfac
    h = nip_spm_solvers(Factors{1}(:,i), L.^2, basis, Qe, 'ARD');
    h = h/max(h);
    for j = 1:size(basis,2)
        Q(:,i) = Q(:,i) + h(j)*basis(:,j);
    end
    %     figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
    %     nip_reconstruction3d(model.cortex,J_rec{i},gca);
end


