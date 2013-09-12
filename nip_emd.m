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

sig1 = nip_energy(Sol1);
sig2 = nip_energy(Sol2);

idx1 = find(sig1);
idx2 = find(sig2);

distance= emd_hat_mex_nes(sig1(idx1),sig2(idx2),aff(idx1,idx2));

end