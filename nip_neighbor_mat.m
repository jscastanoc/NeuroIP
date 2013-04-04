function [Laplacian, QG] = nip_neighbor_mat(cortex)
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
%       QG          -> Aproximation of the green function.
% 
% Additional Comments
% The computation of the Green function is based on a script written by:
% Jose David Lopez - Gareth Barnes - Vladimir Litvak 
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 27 Jan 2013

Nd = size(cortex.vertices,1);
s    = 0.6;
A     = triangulation2adjacency(cortex.faces);
GL    = A - spdiags(sum(A,2),0,Nd,Nd);
GL    = GL*s/2;				% Smoother
Qi    = speye(Nd,Nd);
QG    = sparse(Nd,Nd);
for i = 1:8				% Taylor series approximation
    QG = QG + Qi;
    Qi = Qi*GL/i;
end
QG    = QG.*(QG > exp(-10));		% Eliminate small values
QG    = QG*QG;				% Guarantee positive semidefinite matrix
% clear Qi A GL

options.normalize = 1;
G = compute_mesh_gradient(cortex.vertices,cortex.faces,'distance',options);
Laplacian = sparse(G'*G);
Laplacian = Laplacian + 0.7*abs(min(min(Laplacian)))*eye(size(Laplacian,1));
Laplacian = Laplacian/max(max(Laplacian));
