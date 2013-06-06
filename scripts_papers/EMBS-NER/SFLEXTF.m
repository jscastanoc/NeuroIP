function Q = SFLEXTF(y,L, basis)
% Q = SFLEXTF(y,L, basis)
%
% Juan S. Castano
% 25 May 2013

% Time Frequency representation and Parafac decomposition
Factors = TF_PARAFAC(y);
Nfac = size(Factors{1},2);



% SFLEX Preprocessing
A = L.^2*basis;
nbasis = size(basis,2);
Nd = size(L,2);

% Reconstruction with each factor
Q = zeros(Nd,Nfac);
for i = 1:Nfac
    [xx,status]=dalsql1(ones(nbasis,1), A, Factors{1}(:,i), 0.0007);    
    if sum(find(xx))==0;
        Q(:,i) = ones(Nd,1);
    else
        Q(:,i)  = basis*xx;
    end
%     figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
%     nip_reconstruction3d(model.cortex,J_rec{i},gca);
end



%
% figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
% nip_reconstruction3d(model.cortex,nip_lcmv(model.y,model.L)',gca);
