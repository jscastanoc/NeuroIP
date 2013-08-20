function [Laplacian] = nip_neighbor_mat(cortex)
% [Laplacian, QG] = nip_neighbor_mat(cortex)
% Computes the Laplacian matrix and the green function from the graph
% laplacian.
%
% Input:
%       cortex -> struct. Structure containg the cortex.vertices and 
%               cortex.faces of the graph resulting from the tessellation of
%               the cortex surface.
%
% Output:
%       Laplacian   -> NdxNd. Graph Laplacian.
% 
% Additional Comments
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 27 Jan 2013
if ~isfield(cortex, 'vc') && ~isfield(cortex,'tri')
    cortex.vc = cortex.vertices;
    cortex.vertices = [];
    cortex.tri = cortex.faces;
    cortex.faces = [];
end

Nd = size(cortex.vc,1);

options.normalize = 1;
G = compute_mesh_gradient(cortex.vc,cortex.tri,'distance',options);
Laplacian = sparse(G'*G);
Laplacian = Laplacian + 0.7*abs(min(min(Laplacian)))*eye(size(Laplacian,1));
Laplacian = Laplacian/max(max(Laplacian));
