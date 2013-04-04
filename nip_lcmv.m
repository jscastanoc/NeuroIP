function  Q = nip_lcmv(y,L)
% Q = nip_lcmv(y,L)
%
% Juan S. Castano 
% 12 Mar 2013
Nd = size(L,2);

for i = 1:Nd
    delta(i) = 1/(L(:,i)'*L(:,i));
end
InvCov = spm_inv(y*y');            
for i = 1:Nd
        Q(i) = 1/(L(:,i)'*InvCov*L(:,i));
        Q(i) = Q(i)/delta(i);
end

end