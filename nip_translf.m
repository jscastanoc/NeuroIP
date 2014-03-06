function Ltrans = nip_translf(L)
% function Ltrans = nip_translf(L)
% This function takes a Leadfield matrix and turns it into a lead field
% tensor, where each slice corresponds to the lead field matrix in
% one direction (x y or z). (Works the other way around too!)
% Juan S. Castano C. 
% jscastanoc@gmail.com
% 15 Aug 2013


if ndims(L) == 3
    idx = (1:3:size(L,2));
    Ltrans = reshape(permute(L,[1 3 2]), size(L,1), size(L,2)*3);
elseif ndims(L) == 2
    idx = (1:3:size(L,2)-2);
    Ltrans = zeros([size(L,1) size(L,2)/3 3]);
    for i = 1:3
        Ltrans(:,:,i) = L(:,idx+i-1);
    end
else
    error('L should have at least 2 dimensions')
end
