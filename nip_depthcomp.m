function Lcomp = nip_depthcomp(L,gamma)
% Lcomp = nip_depthcomp(L,gamma)
% Compensate the depth bias on the reconstruction by normalizing the lead
% field matrix. 
%  Input:
%       L       -> Ncx3Nd. Lead field matrix.
%       gamma   -> Scalar. Normalization parameter. 0 means no
%           normalization. 1 means maximum normalization.
%  Output:
%       Lcomp   -> Ncx3Nd. Normalized lead field matrix.
% Juan S. Castano C.
% jscastanoc@gmail.com
% 14 Aug 2013

[Nc Nd] = size(L);

if nargin == 1
    gamma = 0.2;
end

%Normalization of the Lead field matrix to Compensate source depth.
index = (1:3:Nd);
for i = index
    norm_term = sqrt((norm(L(:,i),2)^2+ ...
        norm(L(:,i+1),2)^2 + norm(L(:,i+2),2)^2)^gamma);
    Lcomp(:,i:i+2) = L(:,i:i+2)/norm_term;
end