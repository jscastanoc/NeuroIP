function h = nip_show_3D(cortex, data, a_handle)
% h = nip_show_3D(cortex, data, a_handle)
% Shows the activity 'data' in the volume described by 'mesh'.
% Input:
%       cortex -> struct. Describes the volume to be drawn. Should contain
%               the fields 'faces' and 'vertices' corresponding to the graph 
%               of the tessellated brain surface.
%       data -> Ndx1. Vector with an unnormalized representation of the
%               colors of each vertex mesh.vertices.
%       a_handle -> axes handle. This is where the volume is going to be
%               drawn.
% Output:
%       h -> patch handle.
%
% Additional comments :
% Based on a script written by Eduardo Giraldo.
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 26 Jan 2013


% Map the information contained in data to the corresponding colormap
% codification
data = data-min(data);

x_tik = data;

insig_idx =  find(abs(data) < max(abs(data))*0.1);

nc=256;
[n,T]=size(x_tik);
[v1,vi]=max(x_tik,[],2);
[v2,vi2]=max(v1);
vc= autumn(nc); % If you want to change the colormap used, here is where you should do it
[c,x1] = hist(x_tik(:,vi(vi2)),nc);
di=mean(diff(x1));
vca=zeros(n,3);
k=1;
for k1=1:nc;
    In=find((x_tik(:,k)>=x1(k1)-di) & (x_tik(:,k)<x1(k1)+di));
    vca(In,1:3)=(vc(k1,:)'*ones(1,length(In)))';
end

vca(insig_idx,:) = repmat([255 255 255]/275, length(insig_idx),1);

% Draw (?)
axes(a_handle)
cla
% cortex_smooth = smoothpatch(cortex,1,2,1,0.5);
cortex_smooth = cortex;
h = patch('Faces', cortex_smooth.faces, 'Vertices', cortex_smooth.vertices,'FaceVertexCData',vca,'FaceColor','interp');
axis equal;
axis off;
set(h,'edgecolor','none');
set(h,'AmbientStrength',1,'DiffuseStrength',1.0,'SpecularColorReflectance',0.0)
material dull;
camlight headlight; 
lighting phong
end


