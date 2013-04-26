function [X, act_dip] = nip_simulate_activity(cortex, Laplacian, act_dip, act, t, options)
% X = nip_simulate_activity(cortex, Laplacian, act_dip, act, t, options)
% 
% Input:
%       cortex      -> struct. Describes the volume to be drawn. Should contain
%                   the fields 'faces' and 'vertices' corresponding to the graph 
%                   of the tessellated brain surface.
%		Laplacian   -> NdxNd. Spatial Laplacian matrix (information about dipole Neighbors).
%       act_dip     -> Scalar or Nactx3. If scalar, represent the number of
%                   active sources in the brain, randomly selected.
%                   If matrix, is the coordinates of the active dipoles.
%                   The function looks for the closest dipoles to the given
%                   coordinates.
%       act         -> Same rows (or scalar) than act_dip xNt. Time series
%                   of the activity that's going to be simulated in the 
%                   active dipoles
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
%
% Juan S. Castano
% jscastanoc@gmail.com
% 26 Jan 2013
% TODO: No need for laplacian!

Nd = size(cortex.vertices,1);
Laplacian = speye(Nd);
Nt = length(t);

options.null = 0;

if ~isfield(options, 'sample_all')
    options.sample_all = 0;
end

if isscalar(act_dip)    
        if options.sample_all
            act_dip = randsample(Nd,act_dip)
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

if ischar(act)
    switch act
        case 'real' % This is only for a sample frequency of 256 Hz
			warning('Only valid for a sample frequency of 256 Hz')
            a = [1.0628, 0.000143, -0.000286, -0.42857];
            b = [-0.12];
            X = zeros(Nd,Nt);
            rand_signal = zeros(Nd,Nt);
            rand_signal(act_dip,:) = 0.05*randn(length(act_dip),Nt);
            for i = 3:Nt
                X(:,i) = a(1)*X(:,i-1) + b(1)*Laplacian*X(:,i-1) ...
                    +a(2)*X(:,i-1).^2+a(3)*X(:,i-1).^3+a(4)*X(:,i-2)+rand_signal(:,i);
                if i >Nt/2
                    a(4) = -1;
                    a(1) = 1.3;
                end
            end            
        otherwise
            error('Error: The activity simulator that you selected is not available')
    end
else
    X = zeros(Nd,Nt);
    for i = 1:length(act_dip);
        X(act_dip(i),:) = act(i,:);
    end    
end

end
