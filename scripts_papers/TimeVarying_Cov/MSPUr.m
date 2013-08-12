function Q = MSPUr(yp, L, basis,method)
% Q = MSPUr(yp, L, basis)
%
% Juan S. Castano
% 14 Jun 2013

nbasis = size(basis,2);
Nd = size(L,2);
Nc = size(L,1);
Nr = size(yp,2);
% Reconstruction with each factor
Q = zeros(Nd,Nr);
Qe = eye(Nc);
for i = 1:Nr
    h = nip_spm_priors(yp(:,i), L, basis, Qe,method);
    h = h/max(h);
    for j = 1:size(basis,2)
        Q(:,i) = Q(:,i) + h(j)*basis(:,j);
    end
end