function distance = nip_emd(Sol1, Sol2, aff)
% distance = nip_emd(Sol1, Sol2, aff)
% Computes the earth moving distance between simulations and
% reconstruction.
% Input:
%       Sol1 and Sol2 -> 3Ndx1. activity in each dipole.
%       aff -> NdxNd. Affinity/distance matrix between all the dipoles.
% Output:
%       distance -> Scalar. Earth mover's distance.
%
% Juan S. Castano 
% jscastanoc@gmail.com
% 2 Sep 2013.

Nd = length(Sol1);

data_m = zeros(Nd/3,1);
for i = 1:Nd/3
    data_m(i) = sqrt(sum(Sol1((i-1)*3+1:(i-1)*3+3).^2));
end
sig1 = data_m;

data_m = zeros(Nd/3,1);
for i = 1:Nd/3
    data_m(i) = sqrt(sum(Sol2((i-1)*3+1:(i-1)*3+3).^2));
end
sig2 = data_m;

idx1 = find(sig1);
idx2 = find(sig2);
try
    distance= emd_hat_mex_nes(sig1(idx1),sig2(idx2),aff(idx1,idx2));
catch
    distance= emd_hat_mex(sig1,sig2,aff);
end
end