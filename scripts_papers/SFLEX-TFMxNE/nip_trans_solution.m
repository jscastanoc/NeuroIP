function J_rec = nip_trans_solution(S)
% function J_rec = nip_trans_solution(S)
%
% Juan S. Castanoo C.
% jscastanoc@gmail.com
% 14 Mar 2013
if ndims(S) == 3
    J_rec = permute(S,[3 1 2]);
    J_rec = nip_translf(J_rec)';
elseif ndims(S) == 2
    S = nip_translf(S');
    J_rec = permute(S,[2 1 3]);
else

end