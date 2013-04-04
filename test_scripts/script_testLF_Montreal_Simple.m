close all;clear; clc;

nip_init();
load data/montreal2000_10-10


model.Nd = size(L,2)
model.Nc = size(L,1);
model.Nt = 1;
model.vol = cortex_mesh;

model.L = L;

vol = model.vol;

x = zeros([model.Nd , 3]);
[idx_act , d] = dsearchn(vol.vertices,[200 -50 100]);
x(idx_act,:) = [1 1 1];

point = [ 0 0 0];
dists = dist([vol.vertices; point]' );
% plot(dists(end,:))
[~, idx] = sort(dists(end,:),'descend');

% temp = 1.5;
% x(idx(1:end/temp),:) = ones(size(x(idx(1:end/temp),:)));


% scatter3(model.pos(idx(1:end/temp),1),model.pos(idx(1:end/temp),2),model.pos(idx(1:end/temp),3))
% scatter3(model.pos(:,1),model.pos(:,2),model.pos(:,3))



figure(1)
scatter3(point(1),point(2),point(3))
h = patch('Faces', vol.faces, 'Vertices', vol.vertices,'FaceVertexCData',x,'FaceColor','interp');

set(h,'edgecolor','none')
axis equal
axis off
view([0 90]);

x = x(:,1);
x = zeros([ model.Nd 1]);
x(idx_act) = 1;
% y = model.L*x;
y= L*x;
y = (y-min(y))/max(y-min(y));
% figure
% h = patch('Faces', vol_s.tri, 'Vertices', vol_s.pnt*0.9);
% set(h,'facealpha',0.5);

hold on
% elect_pos(7,:) = [elect_pos(7,1)+23,elect_pos(7,2)+8,elect_pos(7,3)+25];
elect_pos = elec.elecpos;
vol_s.pnt = head{1}.vertices;
vol_s.tri = head{1}.faces;


[idx_act , d] = dsearchn(elect_pos,vol_s.pnt);

for i = 1:length(idx_act);
y_col(i,:) = [y(idx_act(i)),y(idx_act(i)),y(idx_act(i))];
end
%y_col = ones([length(idx_act) 3]);
y = y_col;
y_col = (y-min(min(y)))/max(max(y-min(min(y))));
figure(2)
h = patch('Faces', vol_s.tri, 'Vertices', vol_s.pnt,'FaceVertexCData',y_col,'FaceColor','interp');
set(h,'edgecolor','none')
hold on
scatter3(elect_pos(:,1),elect_pos(:,2),elect_pos(:,3))
axis equal
axis off
view([0 90]);
