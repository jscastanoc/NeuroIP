function  q = nip_lcmv(y,L)
% Q = nip_lcmv(y,L)
% Compute the lcmv (minimum variance) Beamformer 
% Input:
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
%
% Output:
%       q -> Ndx1. Vector containing the beamformer solution
% Juan S. Castano 
% 12 Mar 2013
Nd = size(L,2);

for i = 1:Nd
    delta(i) = 1/(L(:,i)'*L(:,i));
end
InvCov = spm_inv(y*y');      
q = zeros(Nd,1);
for i = 1:Nd
        q(i) = 1/(L(:,i)'*InvCov*L(:,i));
        q(i) = q(i)/delta(i);
end

end