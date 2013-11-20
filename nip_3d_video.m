function nip_3d_video(cortex,data,y,t,varargin)
% function nip_3d_video(cortex,data,y,t,options)
% Shows a video of the activity in contained in data and the time series of
% the electrode measurements.
% Input:
%       cortex -> struct. Describes the volume to be drawn. Should contain
%               the fields 'faces' (or tri) and 'vertices' (or vc) corresponding to the graph 
%               of the tessellated brain surface.
%       data -> 3NdxNt. Time series of the activity in the brain. The format/ordering of the data is
%               (x1,y1,z1,....xNd,yNd,zNd)' for each time instant. where xn yn and zn are the activity
%               of the n-th dipole in each coordinate.
%       y -> NcxNt. Time series of the EEG (electrodes' measurements).
%       t -> 1xNt. Time vector (in secs of ms) with which y was
%               recorded/simulated.
%       options -> struct.
%           	figure. figure id where the video is going to be shown.
% Output:
%
% Juan S. Castano C.
% jscastanoc@gmail.com
% 2 sep 2013


options = varargin{1};
if ~isfield(varargin{1},'figureid'); 
    options.figureid = figure; 
end
if ~isfield(varargin{1},'view'); 
    options.view = [-90,0]; 
end

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

data_m =abs(data_m)- min(abs(data_m(:)));
data_m = data_m/max(data_m(:));

nc=256;
vc= jet(nc); % If you want to change the colormap used, here is where you should do it

figure(options.figureid);
subplot(2,1,1)

if ~isempty(t)&&~isempty(y)
plot(t,y);
hold on
time_mkr = plot([0 0], [max(y(:)) min(y(:))],'r');
end

for i= 1:size(data,2)
    if ~isempty(t)&&~isempty(y)
    set(time_mkr, 'XData', [t(i) t(i)]);
    end
    data = floor(data_m(:,i)*(nc-1))+1;
    vca = vc(data,:);
    
    
%     if i == 1
        subplot(2,1,2)
        axes(gca);
        cla
        h = patch('Faces', cortex.faces, 'Vertices', cortex.vertices,'FaceVertexCData',vca,'FaceColor','interp');
        colormap jet
        caxis([minc maxc])
        axis equal;
        axis off;
        set(h,'edgecolor','none');
        set(h,'AmbientStrength',1,'DiffuseStrength',1.0,'SpecularColorReflectance',0.0)
        material dull;
        lighting phong;
        view(options.view)
%     else
%         set(h,'FaceVertexCData',vca)
%     end
    camlight headlight;
    pause(0.05)
%     i
end
end
