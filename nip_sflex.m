function [J_est, extras] = nip_sflex(y, L, basis, reg_par)
%  [J_est, extras] = nip_sflex(y, L, basis)
% Implements "Large-scale EEG/MEG source localization with spatial
% flexibility." by Haufe et al 2011
%
% Input:
%       y -> NcxNt. Matrix containing the data,
%       L -> NcxNd. Lead Field matrix
%       basis-> NdxNs. Matrix containing the spatial basis functions.
%       reg_par-> Scalar. Regularization parameter (1e-6 by default)
% Output:
%       J_rec -> NdxNt. Reconstructed activity (solution)
%       extras.regpar -> Scalar. Optimum regularization parameter
%
% Additional comments: Uses the DAL optimization toolbox.
% 
% Juan S. Castano C.
% 13 June 2013


[Nc Nt] = size(y);
Nd = size(L,2);

if nargin <=3
    reg_par = 1e-6;
end


% --- Not sure if I should do this here (you'll need a lot of ram if data
% is big enough
A = sparse(kron(speye(Nt), L*basis)); 

nbasis = size(basis,2);

[xx,status]=dalsql1(zeros(nbasis*Nt,1), A, y(:), reg_par);

J_est = basis*reshape(xx,nbasis,Nt);
extras =[];
end