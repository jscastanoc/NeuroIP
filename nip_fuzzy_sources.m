function aff = nip_fuzzy_sources(cortex, sigma, varargin)
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
% Juan S. Castanoo C.
% jscastanoc@gmail.com
% 14 Mar 2013


if ~isfield(cortex, 'vc') && ~isfield(cortex,'tri')
    cortex.vc = cortex.vertices;
    cortex.vertices = [];
    cortex.tri = cortex.faces;
    cortex.faces = [];
end


if isempty(varargin)
    varargin{1} = [];
end
if ~isfield(varargin{1},'dataset')
    varargin{1}.dataset = [];
end
if  ~isfield(varargin{1},'save')
    varargin{1}.save = 0;
end


if ( isempty(varargin{1}.dataset) || ~varargin{1}.save )
    aff = graphrbf(cortex);
else
    filename = strcat('VC',num2str(size(cortex.vc,1)),varargin{1}.dataset,'.mat');
    if exist(filename)
        load(filename);
    else
        aff = graphrbf(cortex);
        if varargin{1}.save
            save(filename, 'aff');
        end
    end
    
end

aff = exp(-aff.^2/sigma^2);