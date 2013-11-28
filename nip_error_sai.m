function [out, Ms, Mr] = nip_error_sai(cortex, J_sim,J_rec, r)
% [sai, Ms, Mr]  = nip_error_sai(cortex, J_sim,J_rec, r)
%
% Computes an spatial accuracy index to evaluate the spatial quality of an
% inverse solution.
%
% Input:
%       cortex -> struct. Describes the volume to be drawn. Should contain
%               the fields 'faces' and 'vertices' corresponding to the graph 
%               of the tessellated brain surface.
%       J_sim -> NdxNt. Simulated brain activity.
%       J_rec -> NdxNt. Reconstructed brain activity.
%       r -> scalar. If a local maximum in the reconstructed activity is
%           with a radius r (using geodesic distance) of a local maximum
%           in the simulated activity, then it is considered a True Positive
% Output:
%       out -> struct. Contains sai, spec (specificity), and sens (sensitivity)
%       Ms ->   Ndx1. Indexes local maxima of the simulated activity 
%                   Contains 1 if the corresponding dipole is a
%                   local maximum. 0 otherwise.
%       Mr ->   Ndx1. Indexes local maxima of the reconstructed activity 
%                   Contains 1 if the corresponding dipole is a
%                   local maximum. 0 otherwise.
%
% Additional Comments: See Belardinelli et al 2012 for further information about
% this error measurement (Source Reconstruction Accuracy of MEG and EEG...)
%
% TO FIX: At the moment, it "tries" to normalize with respect to the
% average distance between neighboring dipoles. But it does not work
% completely good.
%
% Juan S. Castano C.
% 20 May 2013.


Nt = size(J_sim,2);

if ~isfield(cortex, 'vc') && ~isfield(cortex,'tri')
    cortex.vc = cortex.vertices;
    cortex.vertices = [];
    cortex.tri = cortex.faces;
    cortex.faces = [];
end

Nd = size(cortex.vc,1);
file_name = strcat(fileparts(which('nip_init')),'/data/','dist_mat',num2str(Nd),'.mat');


D = nip_fuzzy_sources(cortex,[],struct('dataset','montreal','save',true,'calc','dist'));

GeoD = D;

% Create spatial "masks"
sp_tol = r; % If the energy of a dipole is the biggest within an area of sp_tol, then it is a local maxima
idx = find(D > sp_tol);
D(idx) = 0;
idx = find(D);
D(idx) = 1;
D = D+speye(Nd);

% Average activity and Normalize

Es = nip_energy(mean(J_sim.^2,2));
Es = Es-min(Es);
Es = Es/max(Es);
Er = nip_energy(mean(J_rec.^2,2));
Er = Er-min(Er);
Er = Er/max(Er);

% Threshold values of the average (denoise(?))
thr = 0.02;
idx = find(Es < thr);
Es(idx) = 0;
idx = find(Er < thr);
Er(idx) = 0;

% Find local maxima
Ms = sparse(Nd,1);
Mr = sparse(Nd,1);
for i = find(Es)'  
    % Check if activity in the i-th vertex is a local maxima in the
    % simulated activity
    Cpatch = D(:,i).*Es;
    idx = find(Cpatch > Es(i));
    if sum(idx)==0
        Ms(i) = 1;
    end
end
for i = find(Er)'  
    % The same thing for the reconstructed activity    
    Cpatch = D(:,i).*Er;
    idx = find(Cpatch > Er(i));
    if sum(idx)==0
        Mr(i) = 1;
    end
end


TP = 0; % True positives
FP = 0; % False positives
D = GeoD;
idx = find(D > r);
D(idx) = 0;
idx = find(D);
D(idx) = 1;
D = D+speye(Nd);
for i = find(Mr)'
    if sum(D(:,i).*Ms) > 0
        TP = TP + 1;
    else
        FP = FP + 1;
    end    
end

TN = 0;
FN = 0;
for i = find(~Mr)'
    if sum(D(:,i).*~Ms) > 0
        TN = TN +1;
    else
        FN = FN +1;
    end
end

specificity = TN/(TN+FP);
sensitivity = TP/(TP+FN);

sai = TP/(TP+FP);

out.sai = sai;
out.spec = specificity;
out.sens = sensitivity;
end

