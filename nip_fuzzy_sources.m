function aff = nip_fuzzy_sources(cortex, sigma, file_name)
% aff = nip_fuzzy_sources(cortex, sigma)
% This function returns a matrix that contains information about the
% distance between points in a graph. It places a gaussian with variance = sigma
% in each vertex
% Input:
%       cortex -> Struct. Structure containing the vertices and faces of
%       the graph
%       sigma -> Scalar. Variance of the gaussian placed at each vertex
% Output:
%       aff -> NdxNd. Symmetrical matrix in which the i-th column
%       represents the gaussian placed a round the i-th vertex.
%
% Additional comments: This function uses the graph toolbox to compute the
% distance between each vertex.
%
% Juan S. Castaño C.
% 14 Mar 2013
if isfield(cortex, 'vc') && isfield(cortex,'tri')
    cortex.vertices = cortex.vc;
    cortex.vc = [];
    cortex.faces = cortex.tri;
    cortex.tri = [];
end

Nd = num2str(size(cortex.vertices,1));
A    = sparse(triangulation2adjacency(cortex.faces));

% Search for a file with the precompute geodesic distances. If not found,
% computes them and saves them in a file (This WILL take a while, grab a
% snickers).
if nargin < 3
    file_name = strcat(fileparts(which('nip_init')),'/data/','dist_mat',num2str(Nd),'.mat');
end
if exist(file_name,'file')
    load(file_name)
else
    D   = compute_distance_graph(A);
    save(file_name,'D');
end

% Normalization (so other functions work for meshes of different sizes)
dd = pdist(cortex.vertices,'euclidean');
dd = squareform(dd);
dist_neighbor = mean(mean(sparse(A.*dd)));
nD = dist_neighbor*D;
nD = nD/max(max(nD));

% Affinity /distance matrix with gaussians aroung each vertex
aff = exp(-nD/sigma);