function [J_rec, extras, invT] = nip_loreta(y, L, varargin)
% [J_rec, extras] = nip_loreta(y, L, varargin) 
% Calculate the inverse problem solution using the minimum norm approach
% Input:
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
%       Q -> NdxNd Covariance matrix of the sources (structure, it gets scaled in the algorithm).
% Output:
%       J_rec -> NdxNt. Reconstructed activity (solution)
%       extras.regpar -> Scalar. Optimum regularization parameter
%
% Additional comments: The optimum regularization parameter is calculated
% using the general cross validation function (GCV). See Grech et al. 2008
% for further information.
%
% Juan S. Castano C. 
% jscastanoc@gmail.com
% 26 Jan 2013

Nd = size(L,2)
p = inputParser;
def_Q = speye(Nd);
def_Winv = [];
addParamValue(p,'cov',def_Q);
addParamValue(p,'Winv',def_Winv);
parse(p,varargin{:})
options = p.Results;

warning off % The inv_Lap calculation generates an RCOND problem. Why? don't know. FIX!

% Pre calculation of some constants to speed up the optimization process.
Q = options.cov;
inv_Lap = Q;
eye_Nc = speye(size(L,1));
iLAP_LT = inv_Lap*L';

fprintf('Computing GCV for LORETA... ')
tic
% Get the optimal regularization parameter
gcv_fun = @(alpha) gcv(y,L,alpha, inv_Lap, iLAP_LT, eye_Nc);
% options = optimset('Display','iter','tolX',1e-6);
optionsopt = optimset('tolX',1e-6);
alpha = fminsearch(gcv_fun, 0.5,optionsopt);

% Solution
invT = iLAP_LT/(L*iLAP_LT+abs(alpha)*eye_Nc);
J_rec = invT*y;

if ~isempty(options.Winv)
    J_est = nip_translf(J_rec');
    siJ = size(J_est);
    J_rec = permute(reshape(full(reshape(permute(J_est, [1 3 2]), siJ(1), [])*options.Winv), siJ(1), siJ(3), siJ(2)), [1 3 2]);
    J_rec = nip_translf(J_rec)';
end

extras.regpar = alpha^2;
fprintf('done! \nElapsed time: %.2d secs \n', toc)
end


function gcv_val = gcv(y,L,alpha, inv_Lap, iLAP_LT, eye_Nc)

T = iLAP_LT/(L*iLAP_LT+abs(alpha)*eye_Nc);
x_est = T*y;
A = norm(L*x_est - y,2);
gcv_val = sum(diag(A*A'))/trace((eye_Nc-L*T))^2;

end
