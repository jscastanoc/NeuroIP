function aff = nip_fuzzy_sources(cortex, sigma)
%
% Juan S. Castaño C.
% 14 Mar 2013
Nd = num2str(size(cortex.vertices,1));

file_name = strcat(fileparts(which('nip_init')),'/data/','dist_mat',num2str(Nd),'.mat');
if exist(file_name,'file')
    load(file_name)
else
    A    = triangulation2adjacency(cortex.faces);
    D   = compute_distance_graph(A);
    save(file_name,'D');
end

aff = exp(-D/sigma);