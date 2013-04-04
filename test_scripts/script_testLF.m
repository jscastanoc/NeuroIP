clear; clc;


load data/montreal4000_10-10
nip_init();



% [vol sampled] = nip_subsample_mesh(vol, 4000);
vol = cortex_mesh
cfg.L =L
cfg.t = 0;
model = nip_create_model(cfg);
clear cfg;

x = zeros([model.Nd , 3]);
aux = vol.vertices;
[idx_act , d] = dsearchn(aux,[100 0 0]);
x(idx_act,:) = [1 1 1];

point = mean(vol.vertices);
point(3) = point(3) - point(3)*0.2;
dists = dist([vol.vertices; point]' );
% plot(dists(end,:))
[~, idx] = sort(dists(end,:),'descend');

temp = 1.6;
x(idx(1:end/temp),:) = ones(size(x(idx(1:end/temp),:)));


% scatter3(model.pos(idx(1:end/temp),1),model.pos(idx(1:end/temp),2),model.pos(idx(1:end/temp),3))
% scatter3(model.pos(:,1),model.pos(:,2),model.pos(:,3))



figure(1)
scatter3(point(1),point(2),point(3))
patch('Faces', vol.faces, 'Vertices', vol.vertices,'FaceVertexCData',x,'FaceColor','interp');
axis equal
axis off
view([0 90]);
