function [reduced_vol sampled_indx] = nip_subsample_mesh(vol, obj_vertices)
% [reduced_vol sampled_indx] = nip_subsample_mesh(vol, obj_vertices)
% Subsamples the volume described by vol.vertices and vol.faces. 
% Input:
%       vol     -> Struct. Contains vertices and faces describing the
%               original volume
%       obj_vertices -> Scalar. Number of vertices in the subsampled
%               volume. Note that it is possible to obtain a couple more or less vertices
%               in the resulting volume.
%
% Output:
%       reduced_vol -> Struct. Contains vertices and faces describing the
%       subsampled volume.
%       sampled_indx -> Nobjx1. Contains the indices of the original
%       vertices sampled in the reduced volume
%
% Additional Comments:
% This function is a wrapper of the Graph Toolbox function "reducepatch"
%
% Juan S. Castano C.
% 19 Feb 2013

if obj_vertices > size(vol.vertices,1)
    error('Desired vertices is greater than the current number of vertices')
else
    R = obj_vertices/size(vol.vertices,1);
    reduced_vol = reducepatch(vol,R);
end

if nargout ==2
    sampled_indx = find(sum(ismember(vol.vertices,reduced_vol.vertices),2) >=2);
end

end
