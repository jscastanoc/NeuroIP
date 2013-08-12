function [X, act_dip] = nip_simulate_activity(cortex, act_dip, act, dir, t, options)
% X = nip_simulate_activity(cortex, Laplacian, act_dip, act, t, options)
% 
% Input:
%       cortex      -> struct. Describes the volume to be drawn. Should contain
%                   the fields 'faces' and 'vertices' corresponding to the graph 
%                   of the tessellated brain surface.
%       act_dip     -> Scalar or Nactx3. If scalar, represent the number of
%                   active sources in the brain, randomly selected.
%                   If matrix, is the coordinates of the active dipoles.
%                   The function looks for the closest dipoles to the given
%                   coordinates.
%       act         -> cell of size 3, each cell with same rows than act_dip xNt if act_dip is Nactx3. act_dipxNt if act_dip is an integer. Time series
%                   of the activity that's going to be simulated in the 
%                   active dipoles
%       dir         -> Nactx3. Direction of the dipole activity (vector).
%       t           -> Ntx1. Row vector with the time vector.
%       options     -> struct. Struct containing several options for the
%                   simulation of the activity.
%                   options.sample_all: The population to sample randomly
%                   placed sources is the whole brain, or just the most
%                   superficial sources.
%       
% Output: 
%       X           -> NdxNt. Time series of all the dipoles, This can be used
%                   as ground truth when evaluating the inversion
%                   algorithms.
%       act_dip     -> Nactx1. Index(es) of the vertices corresponding to
%                   active dipoles
% Juan S. Castano
% jscastanoc@gmail.com
% 26 Jan 2013

if isfield(cortex, 'vc') && isfield(cortex,'tri')
    cortex.vertices = cortex.vc;
    cortex.faces = cortex.tri;
end

Nd = size(cortex.vertices,1);
Nt = length(t);

options.null = 0;

if ~isfield(options, 'sample_all')
    options.sample_all = 0;
end

% Normalize orientation vector
norm_dir = sqrt(sum(dir.^2,2));
norm_dir = repmat(norm_dir,1,3);
dir = dir./norm_dir;


if isscalar(act_dip)    
        if options.sample_all
            act_dip = randsample(Nd,act_dip);
        else   % CAUTION!! this was only design for the montreal database
                % You pick a point that's in the center of the brain and
                % sampled only the dipoles that are farther from that point
            point = mean(cortex.vertices);
            point(3) = point(3) - point(3)*0.2;
            dists = dist([cortex.vertices; point]' );
            [~, idx] = sort(dists(end,:),'descend');
            temp = 1.6;
            act_dip = randsample(Nd/temp,act_dip);
            act_dip = idx(act_dip);
        end
elseif ismatrix(act_dip)
    dists = dist([cortex.vertices; act_dip]' );
    nact = size(act_dip,1);
    act_dip = [];
    for i = 1:nact
        [~, act_dip(i)] = min(dists(end-i+1,1:end-nact));
    end        
else
    error('Error: act_dip should be a Nx3 matrix or a scalar')
end

X = zeros(Nd*3,Nt);
for i = 1:length(act_dip);
    for j = 1:3
        X((act_dip(i)-1)*3+j,:) = act(i,:)*dir(j);
    end
end


end
