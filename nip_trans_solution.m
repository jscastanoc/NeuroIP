function J_rec = nip_trans_solution(S)
% function J_rec = nip_trans_solution(S)
% This function takes the solution as produced by the inverse solution of
% the toolbox used at TU Berlin and returns it a format compatible with
% this toolbox. Works the other way around too
% TU Berlin format:
%       z1 
%      /z2 
%    y1/. 
%   /y2 .
% x1/.
% x2 .
% .
% .
%
% This toolbox:
% x1
% y1
% z1
% x2
% y2
% z2 
% .
% .
% Juan S. Castanoo C.
% jscastanoc@gmail.com
% 14 Mar 2013
[Nd,Nt,~] = size(S);
if ndims(S) == 3
    J_rec = permute(S,[2 1 3]);
    J_rec = nip_translf(J_rec)';
elseif ndims(S) == 2
    S = nip_translf(S');
    J_rec = permute(S,[2 1 3]);
else

end