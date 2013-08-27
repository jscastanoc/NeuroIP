function nip_3d_video(cortex,data,y,t,c_fig)

J_rec = data;
% J_rec = abs(J_rec-min(J_rec(:)));
% J_rec = J_rec/(max(J_rec(:)));

if isfield(cortex, 'vc') && isfield(cortex,'tri')
    cortex.vertices = cortex.vc;
    cortex.faces = cortex.tri;
end

[Nd, Nt] = size(data);
if isempty(t)
    t = 1:1:Nt;
end

% Compute the magnitude of the activity in each dipole
data_m = zeros(Nd/3,Nt);
for j = 1: length(t)
    for i = 1:Nd/3
        data_m(i,j) = sqrt(sum(data((i-1)*3+1:(i-1)*3+3,j).^2));
    end
end


minc= min(data_m(:));
maxc = max(data_m(:));

% data_m = J_rec - repmat(mean(J_rec,2),1,Nt);
% data_m = J_rec - min(J_rec(:));
% data_m = data_m/max(abs(data_m(:)));

data_m =abs(data_m)- min(abs(data_m(:)));
data_m = data_m/max(data_m(:));

nc=256;
vc= jet(nc); % If you want to change the colormap used, here is where you should do it

figure(c_fig)
subplot(2,1,1)
plot(t,y);
hold on
time_mkr = plot([0 0], [max(y(:)) min(y(:))],'r');

for i= 1:length(t)
    set(time_mkr, 'XData', [t(i) t(i)]);
    
    % sigidx = find(data_m(:,i) > max(data_m(:))*0.05);
    data = floor(data_m(:,i)*(nc-1))+1;
%     vca = (255/275)*ones(Nd/3,3);
    vca = vc(data,:);
    
    % vca(insig_idx,:) = repmat([255 255 255]/275, length(insig_idx),1);
    
    subplot(2,1,2)
    axes(gca);
    cla
    cortex_smooth = cortex;
    h = patch('Faces', cortex_smooth.faces, 'Vertices', cortex_smooth.vertices,'FaceVertexCData',vca,'FaceColor','interp');
    colorbar
    colormap jet    
    caxis([minc maxc])
    
    axis equal;
    axis off;
    set(h,'edgecolor','none');
    set(h,'AmbientStrength',1,'DiffuseStrength',1.0,'SpecularColorReflectance',0.0)
    material dull;
    camlight headlight;
    lighting phong;
    pause(0.01)
end
end