function Q = SFLEXUr(yp, L,basis)
% function Q = SFLEXUr(yp, L,basis)
%
% Juan S. Castano
% 14 Jun 2013

Nr = size(yp,2);
% SFLEX Preprocessing
A = L*basis;
nbasis = size(basis,2);
Nd = size(L,2);

% Reconstruction with each factor
Q = zeros(Nd,Nr);
for i = 1:Nr
    [Q(:,i),~] = nip_sflex(yp(:,i), L, basis);
    if sum(find(Q(:,i)))==0;
        Q(:,i) = ones(Nd,1);
    end
%     figure('Units','Normalized','Position',[0.3 0.1 0.3 0.3])
%     nip_reconstruction3d(model.cortex,J_rec{i},gca);
end
