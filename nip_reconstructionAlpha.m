function h = nip_reconstructionAlpha(cortex, data, varargin)
% h = nip_reconstructionAlpha(cortex, data, a_handle)
% Shows the activity 'data' in the volume described by 'cortex' code as alpha values.
% Input:
%       cortex -> struct. Describes the volume to be drawn. Should contain
%               the fields 'faces' (or tri) and 'vertices' (or vc) corresponding to the graph 
%               of the tessellated brain surface.
%       data -> Ndx1. Vector with an unnormalized representation of the
%               colors of each vertex mesh.vertices.
%       options -> varargin Different options for the plot:
%               .view -> 1x2. Vector containing azimut and elevation of the
%                   camera (same as matlab's view fuction).
%               .colormap -> 
%               .thres -> Scalar. Minimum Percentage of energy to display
%                   e.g. if 0.05 then only the dipoles with atleast 5% energy
%                   of the most active dipole are plotted with the colormap,
%                   the "inactive" dipoles are painted gray.
%               .colorbar -> String. 'on' or 'off', controls the display of
%                   the colorbar.
%               .axes -> Axes where the patch is going to be displayed.
%               
% Output:
%       h -> patch handle.
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 26 Jan 2013

options = varargin{1};
if ~isfield(varargin{1},'view'); 
    options.view = [0 0]; 
end

if ~isfield(varargin{1},'alpha'); 
    options.alpha = []; 
end
    
if ~isfield(varargin{1},'colormap'); 
    options.colormap = jet(256); 
end

if ~isfield(varargin{1},'thres'); 
    options.thres = 0.0;
end

if ~isfield(varargin{1},'colorbar'); 
    options.colorbar = 'off';
end

if ~isfield(varargin{1},'axes'); 
    options.axes = gca;
end

if isfield(cortex, 'vc') && isfield(cortex,'tri')
    cortex.vertices = cortex.vc;
    cortex.faces = cortex.tri;
end

Nd = length(data);

% Compute the magnitude of the activity in each dipole
if Nd == size(cortex.vertices,1)
%     data_m = sqrt(data.^2);
    data_m = data;    
else
data_m = zeros(Nd/3,1);
for i = 1:Nd/3
    data_m(i) = sqrt(sum(data((i-1)*3+1:(i-1)*3+3).^2));
end
end

if ~isfield(varargin{1},'crange')
    options.crange = [min(data_m(:)),max(data_m(:))];
end


[crange] = options.crange;



idx = find(data_m < crange(1));
data_m(idx) = crange(1);
idx = find(data_m > crange(2));
data_m(idx) = crange(2);





% if sum(find(data_m <0))
%     data_m = 0.5*data_m./max(abs(data_m));
%     data_m = data_m+0.5;
% end

data_m = data_m./max(abs(data_m));
if ~isempty(find(data_m<0))
   data_m = (data_m+1)/2;
end


insig_idx =  find(abs(data_m) < max(abs(data_m))*options.thres);
sig_idx =  find(abs(data_m) >= max(abs(data_m))*options.thres);

vc = colormap(options.colormap);
nc = length(vc);

datac = floor(data_m*(nc-1))+1;
datac(find(datac > (nc-1))) = nc-1;
datac(find(datac < 1)) = 1;

try
    vca = vc(datac,:);
    vca(insig_idx,:) = repmat([255 255 255]/275, length(insig_idx),1);
    noactdip = false;
catch
    vca = repmat([255 255 255]/275, size(cortex.vertices,1),1);
    noactdip = true;
    warning('No active dipoles mapped')
end
axes(options.axes)
cla
alpha_val = (datac-1)/(nc-1) +0.2; alpha_val(find(alpha_val>1)) = 1;
h = patch('Faces', cortex.faces, 'Vertices', cortex.vertices,'FaceVertexCData',[0,0,0],'FaceVertexAlphaData',(datac-1)/(nc-1),'FaceAlpha','interp','FaceColor','interp');
if ~isempty(options.alpha)
    set(h,'FaceAlpha',options.alpha);
end
colorbar(options.colorbar)
colormap(options.colormap)
view(options.view);

if ~noactdip
    caxis([crange(1)-1e-5 crange(2)+1e-5]) 
end
axis equal;
axis off;
set(h,'edgecolor','none');
set(h,'AmbientStrength',1,'DiffuseStrength',1.0,'SpecularColorReflectance',0.0)
material dull;
camlight headlight; 
lighting phong;
end
