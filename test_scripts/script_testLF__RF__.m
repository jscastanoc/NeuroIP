clear; clc;


load data/default_data
nip_init();



[vol sampled] = nip_subsample_mesh(vol, 4000);
cfg.L = L(:,sampled);
cfg.t = 0;
model = nip_create_model(cfg);
clear cfg;

x = zeros([model.Nd , 3]);
aux = dip_pos(sampled,:);
[idx_act , d] = dsearchn(aux,[100 0 0]);
x(idx_act,:) = [1 1 1];

point = mean(vol.vertices);
point(3) = point(3) - point(3)*1;
dists = dist([vol.vertices; point]' );
% plot(dists(end,:))
[~, idx] = sort(dists(end,:),'descend');

temp = 1.5;
x(idx(1:end/temp),:) = ones(size(x(idx(1:end/temp),:)));


% scatter3(model.pos(idx(1:end/temp),1),model.pos(idx(1:end/temp),2),model.pos(idx(1:end/temp),3))
% scatter3(model.pos(:,1),model.pos(:,2),model.pos(:,3))



figure(1)
scatter3(point(1),point(2),point(3))
patch('Faces', vol.faces, 'Vertices', vol.vertices,'FaceVertexCData',x,'FaceColor','interp');
axis equal
axis off
view([0 90]);

% x = x(:,1);
% y = model.L*x;
% y = (y-min(y))/max(y-min(y));
% % figure
% % h = patch('Faces', vol_s.tri, 'Vertices', vol_s.pnt*0.9);
% % set(h,'facealpha',0.5);
% 
% hold on
% elect_pos(7,:) = [elect_pos(7,1)+23,elect_pos(7,2)+8,elect_pos(7,3)+25];
% % scatter3(elect_pos(:,1),elect_pos(:,2)+7,elect_pos(:,3)-20)
% 
% [idx_act , d] = dsearchn(elect_pos,vol_s.pnt);
% 
% for i = 1:length(idx_act);
% y_col(i,:) = [y(idx_act(i)),y(idx_act(i)),y(idx_act(i))];
% end
% % y_col = ones([length(idx_act) 3]);
% figure(2)
% h = patch('Faces', vol_s.tri, 'Vertices', vol_s.pnt*0.9,'FaceVertexCData',y_col,'FaceColor','interp');
% axis equal
% axis off
% view([0 90]);
