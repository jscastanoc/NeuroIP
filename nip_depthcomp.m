function Lcomp = nip_depthcomp(L,gamma)

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